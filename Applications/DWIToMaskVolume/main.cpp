// VTK includes
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkVersion.h>
#include <vtkImageCast.h>
#include <vtkImageData.h>

// vtkTeem includes
#include <vtkTeemNRRDReader.h>
#include <vtkTeemNRRDWriter.h>

// vtkDMRI includes
#include <vtkTeemEstimateDiffusionTensor.h>
#include <vtkTensorMask.h>

// ITK includes
#include <itkFloatingPointExceptions.h>

// vnl includes
#include <vnl/vnl_double_3.h>


// vtkTeem includes
#include <vtkTeemNRRDReader.h>
#include <vtkTeemNRRDWriter.h>

#define MAT_TOL 1e-6

const std::string dwmri_grad_tag("DWMRI_gradient_");
const std::string dwmri_bvalue_tag("DWMRI_b-value");

//----------------------------------------------------------------------------
int ParseDiffusionInformation(vtkTeemNRRDReader* reader,vtkDoubleArray* gradients_array,vtkDoubleArray* bvalues_array);

bool parse_gradient_key(std::string key, size_t& grad_number, size_t& gradkey_pad_width, std::string err);

bool transformsNotEqual(const vtkMatrix4x4* mat1, const vtkMatrix4x4* mat2);
//--------------------------------------------------------------------------


int ParseDiffusionInformation(vtkTeemNRRDReader* reader,vtkDoubleArray* gradients_array,vtkDoubleArray* bvalues_array) 

{
    // Validate modality tag
    std::string modality(reader->GetHeaderValue("modality"));
    if (modality != "DWMRI")
    {
        std::cout << ("NRRD header missing 'modality: DWMRI' tag!");
        return 0;
    }

    std::map<std::string, std::string> nrrd_keys = reader->GetHeaderKeysMap();

    /*
        Step 1: get DWMRI_b-value
    */
    std::string ref_bvalue_str(reader->GetHeaderValue(dwmri_bvalue_tag.c_str()));
    if (ref_bvalue_str.empty())
    {
        std::cout<< ("Missing 'DWMRI_b-value' tag!");
        return 0;
    }
    double ref_bvalue = atof(ref_bvalue_str.c_str());


    /*
      Step 2: loop over all keys
        - for all DWMRI_gradient_ keys, validate
          - consecutive
          - consistent padding
        - save each gradient to tmp_grads
        - record maximum gradient length
    */
    vtkNew<vtkDoubleArray> tmp_grads;
    tmp_grads->SetNumberOfComponents(3);
    size_t grad_idx = 0;
    size_t gradkey_pad_width = 0;
    double max_grad_norm = 0;
    std::string err;

    std::map<std::string, std::string>::iterator nrrd_keys_iter = nrrd_keys.begin();
    for (; nrrd_keys_iter != nrrd_keys.end(); nrrd_keys_iter++)
    {
        std::string key = nrrd_keys_iter->first;

        if (!parse_gradient_key(key, grad_idx, gradkey_pad_width, err))
        {
            if (err.empty())
            {
                continue;
            }
            else
            {
                std::cout<< err;
                return 0;
            }
        }
        // parse the gradient vector into double[3]
        vnl_double_3 cur_grad(0, 0, 0);
        std::stringstream grad_value_stream(nrrd_keys_iter->second);
        grad_value_stream >> cur_grad[0] >> cur_grad[1] >> cur_grad[2];

        max_grad_norm = std::max(cur_grad.two_norm(), max_grad_norm);
        tmp_grads->InsertNextTuple(cur_grad.data_block());
    }

    assert(grad_idx == (size_t)tmp_grads->GetNumberOfTuples());

    /*
      Step 3: loop over gradients
        - calculate each b-value based on NA-MIC DWI gradient length-encoding
        - then normalize each gradient to unit-length
    */
    bvalues_array->SetNumberOfTuples(tmp_grads->GetNumberOfTuples());
    // calculate the b-values
    for (int i = 0; i < tmp_grads->GetNumberOfTuples(); i++)
    {
        vnl_double_3 cur_grad(0, 0, 0);
        cur_grad.copy_in(tmp_grads->GetTuple3(i));

        // note: this is norm^2, per the NA-MIC NRRD DWI convention
        // http://wiki.na-mic.org/Wiki/index.php/NAMIC_Wiki:DTI:Nrrd_format
        double cur_bval = ref_bvalue * pow(cur_grad.two_norm() / max_grad_norm, 2);
        bvalues_array->SetValue(i, cur_bval);

        // normalize gradient vector to unit-length
        //   must be done *after* bvalue extraction
        cur_grad.normalize();
        tmp_grads->InsertTuple(i, cur_grad.data_block());
    }

    // Step 4: copy tmp_grads to output
    gradients_array->DeepCopy(tmp_grads.GetPointer());
    return 1;
}


bool parse_gradient_key(std::string key, size_t& grad_number, size_t& gradkey_pad_width, std::string err)
{
	std::stringstream err_stream;

	// note: Slicer no longer supports NRRDs with DWMRI_NEX_ keys. This was removed
	//       from Dicom2Nrrd in the following commit:
	//           Slicer3/4 SVN: r26101, git-svn: 63a18f7d6900a
	//       and never un-commented in Dicom2Nrrd (or DWIConvert).
	//       If such a key is found, we print an error and fail.
	if (key.find("DWMRI_NEX_") != std::string::npos)
	{
		err_stream << "DWMRI_NEX_ NRRD tag is no longer supported (since SVN r6101)."
			<< " Please adjust header manually to encode repeat excitations"
			<< " as unique DWMRI_gradient_###N keys, and re-load.";
		err = err_stream.str();
		return false;
	}

	if (key.find(dwmri_grad_tag) == std::string::npos)
	{
		return false;
	}
	// below here key is: DWMRI_gradient_####

	// padding is the extra zeros to give a specific digit count
	//   0001
	//   ^^^  <- zeros here are padding to 4 digits
	// we enforce the constraint that the padding of the grad keys must be consistent
	if (gradkey_pad_width == 0) {
		gradkey_pad_width = key.size() - dwmri_grad_tag.size();
	}
	else if (gradkey_pad_width != key.size() - dwmri_grad_tag.size())
	{
		err_stream << "DWMRI NRRD gradient key-numbers must have consistent padding (####N)"
			<< " Found tag: '" << key << "' but previous key had padding: " << gradkey_pad_width;
		err = err_stream.str();
		return false;
	}

	// slices key string from `dwmri_grad_tag.size()` to `key.size()`.
	//   e.g.: "DWMRI_gradient_000001" -> "000001"
	std::string cur_grad_num_str = key.substr(dwmri_grad_tag.size(), key.size());
	size_t cur_grad_num = atol(cur_grad_num_str.c_str());

	// enforce monotonic order
	if (cur_grad_num != grad_number)
	{
		err_stream << "DWMRI NRRD gradient key-numbers must be consecutive."
			<< " Found tag: '" << key << "' but previous key was: " << grad_number;
		err = err_stream.str();
		return false;
	}

	grad_number += 1;
	return true;
}




bool transformsNotEqual(const vtkMatrix4x4* mat1, const vtkMatrix4x4* mat2)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (fabs(mat1->GetElement(i, j) - mat2->GetElement(i, j)) > MAT_TOL)
            {
                return true;
            }
        }
    }
    return false;
}


int main(int argc, char* argv[])
{
    const std::string inputVolume = "C:/Users/Bxd/Desktop/testhead/dwi.nrrd";
    const std::string inputMaskVolume = "C:/Users/Bxd/Desktop/testhead/mask.nrrd";
    const std::string outputBaseline = "C:/Users/Bxd/Desktop/testhead/baseline.nrrd";
    const std::string outputTensor = "C:/Users/Bxd/Desktop/testhead/dti.nrrd";

    const std::string estimationMethod = "LS";//"WLS"
    bool ShiftNegativeEigenvalues = false;
    bool applyMask = false;



    vtkNew<vtkTeemNRRDReader> reader;
    reader->SetFileName(inputVolume.c_str());
    reader->Update();
    if (reader->GetReadStatus())
    {
        std::cerr << argv[0] << ": Error reading Diffusion file" << std::endl;
        return EXIT_FAILURE;
    }

    vtkNew<vtkDoubleArray> bValues;
    vtkNew<vtkDoubleArray> grads;

    if (!ParseDiffusionInformation(reader.GetPointer(), grads.GetPointer(), bValues.GetPointer()))
    {
        std::cerr << argv[0] << ": Error parsing Diffusion information" << std::endl;
        return EXIT_FAILURE;
    }

    vtkNew<vtkTeemEstimateDiffusionTensor> estim;

    estim->SetInputConnection(reader->GetOutputPort());
    estim->SetNumberOfGradients(grads->GetNumberOfTuples());
    estim->SetDiffusionGradients(grads.GetPointer());
    estim->SetBValues(bValues.GetPointer());
    estim->SetShiftNegativeEigenvalues(ShiftNegativeEigenvalues);

    // Compute Transformation that brings the gradients to ijk
    // double *sp = reader->GetOutput()->GetSpacing();
    vtkNew<vtkMatrix4x4> mf;
    mf->DeepCopy(reader->GetMeasurementFrameMatrix());
    vtkNew<vtkMatrix4x4> rasToIjkRotation;
    rasToIjkRotation->DeepCopy(reader->GetRasToIjkMatrix());
    // Set Translation to zero
    for (int i = 0; i < 3; i++)
    {
        rasToIjkRotation->SetElement(i, 3, 0);
    }
    // Remove scaling in rasToIjk to make a real rotation matrix
    double col[3];
    for (int jjj = 0; jjj < 3; jjj++)
    {
        for (int iii = 0; iii < 3; iii++)
        {
            col[iii] = rasToIjkRotation->GetElement(iii, jjj);
        }
        vtkMath::Normalize(col);
        for (int iii = 0; iii < 3; iii++)
        {
            rasToIjkRotation->SetElement(iii, jjj, col[iii]);
        }
    }

    vtkNew<vtkTransform> trans;
    trans->PostMultiply();
    trans->SetMatrix(mf.GetPointer());
    trans->Concatenate(rasToIjkRotation.GetPointer());
    trans->Update();

    estim->SetTransform(trans.GetPointer());
    if (estimationMethod == std::string("LS"))
    {
        estim->SetEstimationMethodToLLS();
    }
    else if (estimationMethod == std::string("WLS"))
    {
        estim->SetEstimationMethodToWLS();
    }
    estim->Update();
    vtkImageData* tensorImage = estim->GetOutput();
    // Read the tensor mask
    vtkNew<vtkImageData> mask;
    if (strlen(inputMaskVolume.c_str()) > 0)
    {
        vtkNew<vtkTeemNRRDReader> maskReader;
        maskReader->SetFileName(inputMaskVolume.c_str());
        maskReader->Update();
        if (maskReader->GetReadStatus())
        {
            std::cerr << argv[0] << ": Error reading Mask file" << std::endl;
            return EXIT_FAILURE;
        }

        // Check if the transforms are equal
        if (transformsNotEqual(maskReader->GetRasToIjkMatrix(), reader->GetRasToIjkMatrix()))
        {
            std::cerr << argv[0] << ": Error reading Mask file, wrong coordinate space" << std::endl;
            return EXIT_FAILURE;
        }

        vtkNew<vtkImageCast> cast;
        cast->SetInputConnection(maskReader->GetOutputPort());
        cast->SetOutputScalarTypeToUnsignedChar();
        cast->Update();

        mask->DeepCopy(cast->GetOutput());
        applyMask = true;
    }
    else
    {
        applyMask = false;
    }

    // Mask tensor
    vtkNew<vtkTensorMask> tensorMask;
    tensorMask->SetNumberOfThreads(1);
    if (applyMask)
    {
        tensorMask->SetMaskAlpha(0.0);
        tensorMask->SetInputConnection(estim->GetOutputPort());
        tensorMask->SetMaskInputData(mask.GetPointer());
        tensorMask->Update();
        tensorImage = tensorMask->GetOutput();
    }
    /**/
    // Compute IjkToRas (used by Writer)
    vtkMatrix4x4* ijkToRasMatrix = reader->GetRasToIjkMatrix();
    ijkToRasMatrix->Invert();

    // Don't save the scalars array, only the tensor array.
    // Save tensor
    vtkNew<vtkTeemNRRDWriter> writer;
    tensorImage->GetPointData()->SetScalars(NULL);
    writer->SetInputData(tensorImage);
    writer->SetFileName(outputTensor.c_str());
    writer->UseCompressionOn();
    writer->SetIJKToRASMatrix(ijkToRasMatrix);
    // Compute measurement frame: Take into account that we have transformed
    // the gradients so tensor components are defined in ijk.
    rasToIjkRotation->Invert();
    writer->SetMeasurementFrameMatrix(rasToIjkRotation.GetPointer());
    writer->Write();

    // Save baseline
    vtkNew<vtkTeemNRRDWriter> writer2;
    writer2->SetInputData(estim->GetBaseline());
    writer2->SetFileName(outputBaseline.c_str());
    writer2->UseCompressionOn();
    writer2->SetIJKToRASMatrix(ijkToRasMatrix);
    writer2->Write();
	
	return EXIT_SUCCESS;
}
