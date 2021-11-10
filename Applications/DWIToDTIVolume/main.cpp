// VTK includes
#include <vtkNew.h>
#include <vtkVersion.h>
#include <vtkPoints.h>
#include <vtkImageCast.h>
#include <vtkImageExtractComponents.h>
#include <vtkImageMedian3D.h>
#include <vtkImageSeedConnectivity.h>
#include <vtkImageThresholdConnectivity.h>
#include <vtkImageWeightedSum.h>


// ITK includes
#include <itkFloatingPointExceptions.h>
#include <itkBinaryFillholeImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkMedianImageFilter.h>
#include <itkOtsuMultipleThresholdsImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkSliceBySliceImageFilter.h>

// ITKVtkGlue
#include "itkVTKImageToImageFilter.h"
#include "itkImageToVTKImageFilter.h"

// vnl includes
#include <vnl/vnl_double_3.h>


// vtkTeem includes
#include <vtkTeemNRRDReader.h>
#include <vtkTeemNRRDWriter.h>


const std::string dwmri_grad_tag("DWMRI_gradient_");
const std::string dwmri_bvalue_tag("DWMRI_b-value");

//----------------------------------------------------------------------------
int ParseDiffusionInformation(vtkTeemNRRDReader* reader,vtkDoubleArray* gradients_array,vtkDoubleArray* bvalues_array);

bool parse_gradient_key(std::string key, size_t& grad_number, size_t& gradkey_pad_width, std::string err);
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
int main(int argc, char* argv[])
{

    const std::string inputVolume = "C:/Users/Bxd/Desktop/testhead/dwi.nrrd";
    const std::string outputBaseline = "C:/Users/Bxd/Desktop/testhead/baseline.nrrd";
    const std::string thresholdMask = "C:/Users/Bxd/Desktop/testhead/mask.nrrd";
    const int baselineBValueThreshold = 100;
    bool removeIslands = true;


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

    // Compute the mean baseline image
    vtkNew<vtkImageWeightedSum> imageWeightedSum;
    imageWeightedSum->NormalizeByWeightOn();

    int b0_count = 0;
    for (int bval_n = 0; bval_n < bValues->GetNumberOfTuples(); bval_n++)
    {
        double bvalue = bValues->GetTuple1(bval_n);
        if (bvalue <= baselineBValueThreshold)
        {
            vtkNew<vtkImageExtractComponents> extractComponents;
            extractComponents->SetInputConnection(reader->GetOutputPort());
            extractComponents->SetComponents(bval_n);
            extractComponents->Update();

            imageWeightedSum->AddInputConnection(extractComponents->GetOutputPort());
            imageWeightedSum->SetWeight(b0_count++, 1.);
        }
    }
    imageWeightedSum->Update();

    if (b0_count == 0)
    {
        std::cerr << argv[0] << ": Error parsing Diffusion information, no B0 images" << std::endl;
        return EXIT_FAILURE;
    }

    vtkNew<vtkImageCast> inputCast;
    inputCast->SetOutputScalarTypeToShort();
    inputCast->SetInputConnection(imageWeightedSum->GetOutputPort());
    inputCast->Update();

    typedef itk::Image<short, 3> ImageType;
    typedef itk::VTKImageToImageFilter<ImageType> VTKImageToImageType;

    VTKImageToImageType::Pointer vtkToITK = VTKImageToImageType::New();
    vtkToITK->SetInput(inputCast->GetOutput());

    typedef itk::MedianImageFilter<ImageType, ImageType> medianType;
    medianType::Pointer median = medianType::New();
    median->SetInput(vtkToITK->GetOutput());
    median->Update();

    typedef itk::OtsuMultipleThresholdsImageFilter<ImageType, ImageType> otsuType;
    otsuType::Pointer otsu = otsuType::New();
    otsu->SetInput(median->GetOutput());
    otsu->SetNumberOfThresholds(2);
    otsu->SetNumberOfHistogramBins(128);
    // To request old ITK4 bin behavior, if needed in future.
    //otsu->ReturnBinMidpointOn();
    otsu->Update();

    typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ccType;
    ccType::Pointer cc = ccType::New();
    cc->SetInput(otsu->GetOutput());
    cc->Update();

    typedef itk::RelabelComponentImageFilter<ImageType, ImageType> rlType;
    rlType::Pointer rl = rlType::New();
    rl->SetMinimumObjectSize(10000);
    rl->SetInput(cc->GetOutput());
    rl->Update();

    typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> threshType;
    threshType::Pointer thresh = threshType::New();
    // ITK5 RelabelComponentImageFilter (above) sets background to 1. This is
    // unexpected behavior since it should ignore background (0) voxels.
    // FIX: Check for presence of two output objects.
    // Object with label 2 is second largest (the brain)
    if (rl->GetNumberOfObjects() == 2)
    {
        thresh->SetLowerThreshold(2);
        thresh->SetUpperThreshold(2);

    }
    else
    {
        // ITK4 used this code assuming one output relabeled object
        thresh->SetLowerThreshold(1);
        thresh->SetUpperThreshold(1);
    }
    thresh->SetOutsideValue(0);
    thresh->SetInsideValue(1);
    thresh->SetInput(rl->GetOutput());
    thresh->Update();

    typedef itk::Image<short, 2> BinImageType;
    typedef itk::BinaryFillholeImageFilter<BinImageType> binfillType;
    binfillType::Pointer binFiller = binfillType::New();
    binFiller->SetForegroundValue(1);

    typedef itk::SliceBySliceImageFilter<ImageType, ImageType> sbsType;
    sbsType::Pointer sbsFilter = sbsType::New();
    sbsFilter->SetFilter(binFiller);
    sbsFilter->SetInput(thresh->GetOutput());
    sbsFilter->Update();

    typedef itk::ImageToVTKImageFilter<ImageType> ITKOutType;
    ITKOutType::Pointer itkToVTK = ITKOutType::New();

    vtkNew<vtkImageCast> cast2;
    cast2->SetOutputScalarTypeToUnsignedChar();

    if (removeIslands)
    {
        itkToVTK->SetInput(sbsFilter->GetOutput());
        cast2->SetInputData(itkToVTK->GetOutput());
    }
    else
    {
        itkToVTK->SetInput(thresh->GetOutput());
        cast2->SetInputData(itkToVTK->GetOutput());
    }
    itkToVTK->Update();

    vtkMatrix4x4* ijkToRasMatrix = reader->GetRasToIjkMatrix();
    ijkToRasMatrix->Invert();

    // Save baseline
    vtkNew<vtkTeemNRRDWriter> writer;
    writer->SetInputConnection(imageWeightedSum->GetOutputPort());
    writer->SetFileName(outputBaseline.c_str());
    writer->UseCompressionOn();
    writer->SetIJKToRASMatrix(ijkToRasMatrix);
    writer->Write();

    // Save mask
    vtkNew<vtkTeemNRRDWriter> writer2;
    writer2->SetInputConnection(cast2->GetOutputPort());

    writer2->SetFileName(thresholdMask.c_str());
    writer2->UseCompressionOn();
    writer2->SetIJKToRASMatrix(ijkToRasMatrix);
    writer2->Write();
	
	return EXIT_SUCCESS;
}
