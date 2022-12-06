#include "vtkSeedTracts.h"
#include "vtkTRKReader.h"
#include "vtkTRKWriter.h"

//VTK include
#include <vtkNrrdReader.h>
#include <vtkImageChangeInformation.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkPLYWriter.h>
#include <vtkPointData.h>
#include <vtkNIFTIImageReader.h>
#include <vtkTeemNRRDReader.h>
#include <vtkSetGet.h>
#include <vtkPolyDataWriter.h>
#include <vtkPointSet.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkLine.h>
#include <vtkPolyLine.h>
#include <vtkDemandDrivenPipeline.h>
#include <vtkImageCast.h>
#include <vtkImageThreshold.h>
#include <vtkMatrix3x3.h>

#include <iostream>
#include <vector>
#include <array>

void ConvertBetweenRASAndLPS(vtkPointSet* inputMesh, vtkPointSet* outputMesh);
void ReadNRRD(std::string fileName);
void ReadNIFTI(std::string fileName, vtkImageData* image);
void LPS2RAS(vtkSmartPointer<vtkImageData> input, vtkMatrix4x4* lpsToIjkMatrix);
void seedStreamlineFromPoint(vtkSeedTracts* seed, std::vector<std::array<double, 3>>& seedPoint, double regionSize, double sampleStep);
void CreateTractsForOneSeed(vtkSeedTracts* seed, int thresholdMode,
	double stoppingValue,
	double stoppingCurvature,
	double integrationStepLength,
	double minPathLength,
	double regionSize, double sampleStep,
	int maxNumberOfSeeds,
	int seedSelectedFiducials, double spa[3], vtkMatrix4x4* RASToIJK, vtkImageData* imgData)
{
	// 设置 world to tensorIJK
	vtkNew<vtkMatrix4x4> mat;
	//volumeNode->GetRASToIJKMatrix(mat.GetPointer());

	vtkNew<vtkMatrix4x4> tensorRASToIJK;
	tensorRASToIJK->DeepCopy(RASToIJK);

	vtkNew<vtkTransform> trans;
	trans->Identity();
	trans->PreMultiply();
	trans->SetMatrix(tensorRASToIJK.GetPointer());
	// Trans from IJK to RAS
	trans->Inverse();
	// Take into account spacing to compute Scaled IJK
	trans->Scale(1 / spa[0], 1 / spa[1], 1 / spa[2]);
	trans->Inverse();
	seed->SetWorldToTensorScaledIJK(trans.GetPointer());




	//Set tensor Rotation matrix 
	vtkNew<vtkMatrix4x4> tensorRASToIJKRotation;
	tensorRASToIJKRotation->DeepCopy(tensorRASToIJK.GetPointer());

	//Set Translation to zero
	for (int i = 0; i < 3; i++)
	{
		tensorRASToIJKRotation->SetElement(i, 3, 0);
	}
	//Remove scaling in rasToIjk to make a real roation matrix
	double col[3];
	for (int jjj = 0; jjj < 3; jjj++)
	{
		for (int iii = 0; iii < 3; iii++)
		{
			col[iii] = tensorRASToIJKRotation->GetElement(iii, jjj);
		}
		vtkMath::Normalize(col);
		for (int iii = 0; iii < 3; iii++)
		{
			tensorRASToIJKRotation->SetElement(iii, jjj, col[iii]);
		}
	}
	tensorRASToIJKRotation->Invert();


	seed->SetTensorRotationMatrix(tensorRASToIJKRotation.GetPointer());

	//ROI comes from tensor, IJKToRAS is the same
	// as the tensor
	vtkNew<vtkTransform> trans2;
	trans2->Identity();
	trans2->SetMatrix(tensorRASToIJK.GetPointer());
	trans2->Inverse();
	seed->SetROIToWorld(trans2.GetPointer());

	seed->UseVtkHyperStreamlinePoints();

	//
	vtkNew<vtkHyperStreamlineDTMRI> streamer;
	seed->SetVtkHyperStreamlinePointsSettings(streamer.GetPointer());
	seed->SetMinimumPathLength(minPathLength);

	//
	streamer->SetThresholdModeToFractionalAnisotropy();
	streamer->SetStoppingThreshold(stoppingValue);
	streamer->SetRadiusOfCurvature(stoppingCurvature);
	streamer->SetIntegrationStepLength(integrationStepLength);

	// Temp fix to provide a scalar
	vtkImageData* inputTensorField = vtkImageData::SafeDownCast(
		seed->GetInputTensorFieldConnection()->GetProducer()->GetOutputDataObject(0));

	inputTensorField->GetPointData()->SetScalars(imgData->GetPointData()->GetScalars());

	//std::array<double, 3> points{ -1.160122164827107,-29.925899726775967,60.27479726775956 };
	//std::array<double, 3> points{ -12.873638961089242, -31.32079998433062,	73.2644539595872 };
	std::array<double, 3> points{ -16.11382170952004,-19.05866524362014,44.88020394634446 };
	//std::array<double, 3> points{ -14.1395,-64.8273,64.8273 };

	//std::array<double, 3> points{ 0.39,30.21,59.68 };
	std::vector<std::array<double, 3>> sPoints;
	sPoints.push_back(points);
	seedStreamlineFromPoint(seed, sPoints, regionSize, sampleStep);

	//6. Extract PolyData in RAS
	vtkSmartPointer<vtkPolyData> outFibers = vtkSmartPointer<vtkPolyData>::New();
	seed->TransformStreamlinesToRASAndAppendToPolyData(outFibers.GetPointer());
	seed->UpdateAllHyperStreamlineSettings();



	// We explicitly write the coordinate system into the file header.
	const std::string coordinateSystemTag = "SPACE"; // following NRRD naming convention
	const std::string coordinateSystemStr = "RAS";
	// SPACE=RAS format follows Mimics software's convention, saving extra information into the
	// STL file header: COLOR=rgba,MATERIAL=rgbargbargba
	// (see details at https://en.wikipedia.org/wiki/STL_(file_format)#Color_in_binary_STL)
	const std::string coordinateSytemSpecification = coordinateSystemTag + "=" + coordinateSystemStr;

	//vtkSmartPointer<vtkPointSet> meshToWrite = ;
	//vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	//vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	//std::string fiberFileName = "T:/scene/fiber_ASCII.vtk";
	std::string fiberFileName = "C:/Users/Bxd/Desktop/testhead/fiber_one.vtk";
	writer->SetFileName(fiberFileName.c_str());
#if VTK_MAJOR_VERSION >= 9
	// version 5.1 is not compatible with earlier Slicer versions and most other software
	writer->SetFileVersion(42);
#endif
	writer->SetFileType( VTK_BINARY);
	//writer->SetFileType(VTK_ASCII);

	writer->SetInputData(outFibers.GetPointer());
	std::string header = std::string("3D Slicer output. ") + coordinateSytemSpecification;
	writer->SetHeader(header.c_str());
	writer->Write();
	writer->Update();

}

//----------------------------------------------------------------------------
int CreateTractsForLabelMap(
	vtkSeedTracts* seed,
	int ROIlabel,
	int useIndexSpace,
	double seedSpacing,
	int randomGrid,
	double linearMeasureStart,
	int thresholdMode,
	double stoppingValue,
	double stoppingCurvature,
	double integrationStepLength,
	double minPathLength,
	double maxPathLength, vtkMatrix4x4* RASToIJK, vtkImageData* imgData, vtkMatrix4x4* roiRASToIJK,vtkImageData* labelmap)
{
	if (imgData == nullptr || labelmap == nullptr)
		return 0;

	vtkSmartPointer<vtkAlgorithmOutput> ROIConnection;
	vtkNew<vtkImageCast> imageCast;
	vtkNew<vtkDiffusionTensorMathematics> math;
	vtkNew<vtkImageThreshold> th;
	vtkNew<vtkMatrix4x4> ROIRASToIJK;

	//1.set inout
	if (labelmap)
	{
		// cast roi to short data type
		imageCast->SetOutputScalarTypeToShort();
		imageCast->SetInputData(labelmap);

		//Do scale IJK
		double sp[3];
		labelmap->GetSpacing(sp);
		vtkImageChangeInformation* ici = vtkImageChangeInformation::New();
		ici->SetOutputSpacing(sp);
		imageCast->Update();
		ici->SetInputConnection(imageCast->GetOutputPort());
		ici->Update();
		ROIConnection = ici->GetOutputPort();

		// Set up the matrix that will take points in ROI
		// to RAS space.  Code assumes this is world space
		// since  we have no access to external transforms.
		// This will only work if no transform is applied to
		// ROI and tensor volumes.
		//
		ROIRASToIJK->DeepCopy(roiRASToIJK);
	}

	// 2. Set Up matrices
	// 设置 world to tensorIJK
	vtkNew<vtkMatrix4x4> mat;
	//volumeNode->GetRASToIJKMatrix(mat.GetPointer());
	double spa[3];
	imgData->GetSpacing(spa);

	vtkNew<vtkMatrix4x4> tensorRASToIJK;
	tensorRASToIJK->DeepCopy(RASToIJK);

	vtkNew<vtkTransform> trans;
	trans->Identity();
	trans->PreMultiply();
	trans->SetMatrix(tensorRASToIJK.GetPointer());
	// Trans from IJK to RAS
	trans->Inverse();
	// Take into account spacing to compute Scaled IJK
	trans->Scale(1 / spa[0], 1 / spa[1], 1 / spa[2]);
	trans->Inverse();
	seed->SetWorldToTensorScaledIJK(trans.GetPointer());


	vtkNew<vtkMatrix4x4> tensorRASToIJKRotation;
	tensorRASToIJKRotation->DeepCopy(tensorRASToIJK.GetPointer());
	// Set Translation to zero
	for (int i = 0; i < 3; i++)
	{
		tensorRASToIJKRotation->SetElement(i, 3, 0);
	}
	// Remove scaling in rasToIjk to make a real rotation matrix
	double col[3];
	for (int jjj = 0; jjj < 3; jjj++)
	{
		for (int iii = 0; iii < 3; iii++)
		{
			col[iii] = tensorRASToIJKRotation->GetElement(iii, jjj);
		}
		vtkMath::Normalize(col);
		for (int iii = 0; iii < 3; iii++)
		{
			tensorRASToIJKRotation->SetElement(iii, jjj, col[iii]);
		}
	}
	tensorRASToIJKRotation->Invert();
	seed->SetTensorRotationMatrix(tensorRASToIJKRotation.GetPointer());

	vtkNew<vtkTransform> trans2;
	trans2->Identity();
	trans2->PreMultiply();

	// no longer assume this ROI is in tensor space
	// trans2->SetMatrix(tensorRASToIJK.GetPointer());
	trans2->SetMatrix(ROIRASToIJK.GetPointer());
	trans2->Inverse();
	seed->SetROIToWorld(trans2.GetPointer());


	// PENDING: Do merging with input ROI

	seed->SetInputROIConnection(ROIConnection);
	seed->SetInputROIValue(ROIlabel);
	seed->UseStartingThresholdOn();
	seed->SetStartingThreshold(linearMeasureStart);

	if (useIndexSpace)
	{
		seed->SetIsotropicSeeding(0);
	}
	else
	{
		seed->SetIsotropicSeeding(1);
	}

	seed->SetRandomGrid(randomGrid);

	seed->SetIsotropicSeedingResolution(seedSpacing);
	seed->SetMinimumPathLength(minPathLength);
	seed->UseVtkHyperStreamlinePoints();
	vtkNew<vtkHyperStreamlineDTMRI> streamer;
	seed->SetVtkHyperStreamlinePointsSettings(streamer.GetPointer());

	streamer->SetThresholdModeToFractionalAnisotropy();
	streamer->SetStoppingThreshold(stoppingValue);
	streamer->SetMaximumPropagationDistance(maxPathLength);
	streamer->SetRadiusOfCurvature(stoppingCurvature);
	streamer->SetIntegrationStepLength(integrationStepLength);

	// Temp fix to provide a scalar
// seed->GetInputTensorField()->GetPointData()->SetScalars(math->GetOutput()->GetPointData()->GetScalars());

	// 5. Run the thing
	seed->SeedStreamlinesInROI();


	//6. Extract PolyData in RAS
	vtkSmartPointer<vtkPolyData> outFibers = vtkSmartPointer<vtkPolyData>::New();
	seed->TransformStreamlinesToRASAndAppendToPolyData(outFibers.GetPointer());
	seed->UpdateAllHyperStreamlineSettings();



	// We explicitly write the coordinate system into the file header.
	const std::string coordinateSystemTag = "SPACE"; // following NRRD naming convention
	const std::string coordinateSystemStr = "RAS";
	// SPACE=RAS format follows Mimics software's convention, saving extra information into the
	// STL file header: COLOR=rgba,MATERIAL=rgbargbargba
	// (see details at https://en.wikipedia.org/wiki/STL_(file_format)#Color_in_binary_STL)
	const std::string coordinateSytemSpecification = coordinateSystemTag + "=" + coordinateSystemStr;

	//vtkSmartPointer<vtkPointSet> meshToWrite = ;
	//vtkSmartPointer<vtkPLYWriter> writer = vtkSmartPointer<vtkPLYWriter>::New();

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	//vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	//std::string fiberFileName = "T:/scene/fiber_ASCII.vtk";
	std::string fiberFileName = "C:/Users/Bxd/Desktop/testhead/fiber_label_one.vtk";
	writer->SetFileName(fiberFileName.c_str());
#if VTK_MAJOR_VERSION >= 9
	// version 5.1 is not compatible with earlier Slicer versions and most other software
	writer->SetFileVersion(42);
#endif
	writer->SetFileType(VTK_BINARY);
	//writer->SetFileType(VTK_ASCII);

	writer->SetInputData(outFibers.GetPointer());
	std::string header = std::string("3D Slicer output. ") + coordinateSytemSpecification;
	writer->SetHeader(header.c_str());
	writer->Write();
	writer->Update();

	return 1;

}

void seedStreamlineFromPoint(vtkSeedTracts* seed, std::vector<std::array<double, 3>>& seedPoint, double regionSize, double sampleStep)
{
	//seedPoint 
	//x = -1.1601221648271069,y = -29.925899726775967,z = 60.274797267759560
	/*int numberofPoint = seedPoint.size();*/
	for (const auto& sPoint : seedPoint) {
		for (double x = -regionSize / 2.0; x <= regionSize / 2.0; x += sampleStep)
		{
			for (double y = -regionSize / 2.0; y <= regionSize / 2.0; y += sampleStep)
			{
				for (double z = -regionSize / 2.0; z <= regionSize / 2.0; z += sampleStep)
				{
					float newXYZ[3];
					newXYZ[0] = sPoint[0] + x;
					newXYZ[1] = sPoint[1] + y;
					newXYZ[2] = sPoint[2] + z;
					//float* xyz = transFiducial->TransformFloatPoint(newXYZ);
					//Run the thing
					//seed->SeedStreamlineFromPoint(xyz[0], xyz[1], xyz[2]);
					seed->SeedStreamlineFromPoint(newXYZ[0], newXYZ[1], newXYZ[2]);
				}
			}
		}
	}
}

int main(int argc, char* argv[])
{
#if 1
	//0.create vtkSeedTracts
	vtkNew<vtkSeedTracts> seed;
	//1.read nrrd
	//std::string fileName = "T:/scene/myDTI.nrrd";

	std::string fileName = "C:/Users/Bxd/Desktop/testhead/dti.nrrd";
	//std::string labelFilename = "C:/Users/Bxd/Desktop/testhead/mask.nrrd";
    std::string labelFilename = "C:/Users/Bxd/Desktop/testhead/Segmentation-seg_1-label.nrrd";
	//ReadNRRD(("W:/scene/myDTI.nrrd"));

	vtkNew<vtkTeemNRRDReader> reader;

	reader->SetUseNativeOriginOn();

	reader->SetFileName(fileName.c_str());

	// Check if this is a NRRD file that we can read
	if (!reader->CanReadFile(fileName.c_str()))
	{
		std::cout << ("ReadData: This is not a nrrd file");

		return 0;
	}

	reader->UpdateInformation();

	// Check type

	if (!(reader->GetPointDataType() == vtkDataSetAttributes::TENSORS))
	{
		std::cout << ("ReadData: MRMLVolumeNode does not match file kind");
		return 0;
	}

	reader->Update();
	// set volume attributes
	vtkMatrix4x4* mat = reader->GetRasToIjkMatrix();
	vtkNew<vtkMatrix4x4> mat_inver;
	vtkMatrix4x4::Invert(mat, mat_inver.Get());
	mat_inver->Print(std::cout);
	// set measurement frame
	vtkMatrix4x4* mat2;
	mat2 = reader->GetMeasurementFrameMatrix();


	double sp[3]{ 0 };
	double orgin[3]{ 0 };
	int dim[3]{ 0 };

	vtkImageData* imgData = reader->GetOutput();
	imgData->GetSpacing(sp);
	imgData->GetOrigin(orgin);
	imgData->GetDimensions(dim);
	std::cout << "spacing:(" << sp[0] << "," << sp[1] << "," << sp[2] << ")" << std::endl;
	std::cout << "orgin:(" << orgin[0] << "," << orgin[1] << "," << orgin[2] << ")" << std::endl;
	std::cout << "dimensions:(" << dim[0] << "," << dim[1] << "," << dim[2] << ")" << std::endl;


	vtkNew<vtkImageChangeInformation> ici;
	//1.9531249999999991,1.9531249999999936,3.5100013269750843

	ici->SetOutputSpacing(sp);
	/*double ori[3]{ 0.0,0.0,0.0 };
	ici->SetOutputOrigin(ori);*/
	ici->SetInputConnection(reader->GetOutputPort());
	seed->SetInputTensorFieldConnection(ici->GetOutputPort());



	//2. 设置 seed 需要的各种参数
	int  thresholdMode = 0;
	double stoppingValue = 0.250;
	double stoppingCurvature = 0.699;
	double integrationStepLength = 0.50;
	double minPathLength = 20.0;
	double regionSize = 2.500;
	double sampleStep = 1.00;
	int maxNumberOfSeeds = 100;
	int seedSelectedFiducials = 0;
	int displayMode = 1;

	vtkSmartPointer<vtkImageData> labelmap = vtkSmartPointer<vtkImageData>::New();
	//ReadNIFTI(labelFilename, labelmap);

	//vtkNew<vtkNIFTIImageReader> labelmapReader;
	//labelmapReader->SetFileName(labelFilename.c_str());
	//labelmapReader->Update();
	//labelmap->DeepCopy(labelmapReader->GetOutput());


	//vtkSmartPointer<vtkMatrix4x4> rotateMat = vtkSmartPointer<vtkMatrix4x4>::New();
	//rotateMat->SetElement(0, 0, -1);
	//rotateMat->SetElement(1, 1, -1);

	//vtkMatrix4x4::Multiply4x4(rotateMat, labelmap->GetPhysicalToIndexMatrix(), labelmap->GetPhysicalToIndexMatrix());


	vtkNew<vtkTeemNRRDReader> labelReader;
	labelReader->SetFileName(labelFilename.c_str());
	labelReader->SetUseNativeOriginOn();
	labelReader->UpdateInformation();
	labelReader->Update();
	/*vtkMatrix4x4* mat1 = labelReader->GetRasToIjkMatrix();
	mat1->Print(std::cout);*/
	vtkMatrix4x4* roiRAS2Ijk = labelReader->GetRasToIjkMatrix();
	roiRAS2Ijk->Print(std::cout);
	labelmap->DeepCopy(labelReader->GetOutput());

	//vtkMatrix4x4* lps2ijk = labelmap->GetPhysicalToIndexMatrix();
	////LPS2RAS(labelmap, lps2ijk);
	//vtkMatrix4x4* ijk2lpsMatrix = labelmap->GetIndexToPhysicalMatrix();
	//ijk2lpsMatrix->Print(std::cout);
	//vtkMatrix3x3* dir = labelmap->GetDirectionMatrix();
	//dir->Print(std::cout);

	double sp1[3]{ 0 };
	double orgin1[3]{ 0 };
	int dim1[3]{ 0 };
	labelmap->GetSpacing(sp1);
	labelmap->GetOrigin(orgin1);
	labelmap->GetDimensions(dim1);
	std::cout << "spacing1:(" << sp1[0] << "," << sp1[1] << "," << sp1[2] << ")" << std::endl;
	std::cout << "orgin1:(" << orgin1[0] << "," << orgin1[1] << "," << orgin1[2] << ")" << std::endl;
	std::cout << "dimensions1:(" << dim1[0] << "," << dim1[1] << "," << dim1[2] << ")" << std::endl;
	

	// ITK image direction are in LPS space
		// convert from ijkToRas to ijkToLps
	//vtkMatrix4x4* LPSToRASMatrix = vtkMatrix4x4::New();
	//LPSToRASMatrix->Identity();
	//LPSToRASMatrix->SetElement(0, 0, -1);
	//LPSToRASMatrix->SetElement(1, 1, -1);

	//vtkMatrix4x4* ijkToRASMatrix = vtkMatrix4x4::New();
	//vtkMatrix4x4::Multiply4x4(ijk2lpsMatrix, LPSToRASMatrix, ijkToRASMatrix);

	if (labelmap)
	{
		double range[2];
		range[0] = -1;
		range[1] = 1e10;
		labelmap->GetScalarRange(range);

		//int n = range[1] - range[0];
		//以最大的值为索引值
		int useIndexSpace = 1;
		int ROIlabel = range[1];
		double seedSpacing = 2.0000;
		int randomGrid = 0;
		double linearMeasureStart = 0.29999999;
		int thresholdMode = 0;
		double stoppingValue = 0.25;
		double stoppingCurvature = 0.69999999999999996;
		double integrationStepLength = 0.50000000000000000;
		double minPathLength = 20.000000000000000;
		double maxPathLength = 800.00000000000000;

		CreateTractsForLabelMap(seed, ROIlabel, useIndexSpace, seedSpacing, randomGrid, linearMeasureStart, thresholdMode, stoppingValue, stoppingCurvature,
			integrationStepLength, minPathLength, maxPathLength, mat, ici->GetOutput(), roiRAS2Ijk, labelmap);

	}
	else
	{
		CreateTractsForOneSeed(seed, thresholdMode, stoppingValue, stoppingCurvature, integrationStepLength,
			minPathLength, regionSize, sampleStep, maxNumberOfSeeds, seedSelectedFiducials, sp, mat, ici->GetOutput());
	}

	










	std::cout << "hello visualDMRI" << std::endl;
	return 0;
#endif 

#if 0

	vtkSmartPointer<vtkTRKReader> reader = vtkSmartPointer<vtkTRKReader>::New();
	reader->SetFileName("D:/Test/v3d/zzv3d/AF_left.trk");
	reader->SetFileName("D:/Test/v3d/zzv3d/AF_left_myTest.trk");
	reader->Update();
	


	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	//vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	//std::string fiberFileName = "T:/scene/fiber_ASCII.vtk";
	std::string fiberFileName = "D:/Test/v3d/zzv3d/AF_left.vtk";
	writer->SetFileName(fiberFileName.c_str());
#if VTK_MAJOR_VERSION >= 9
	// version 5.1 is not compatible with earlier Slicer versions and most other software
	writer->SetFileVersion(42);
#endif
	writer->SetFileType(VTK_BINARY);
	//writer->SetFileType(VTK_ASCII);

	writer->SetInputData(reader->GetOutput());
	//std::string header = std::string("3D Slicer output. ") + coordinateSytemSpecification;
	//writer->SetHeader(header.c_str());
	writer->Write();
	writer->Update();


	std::string trkFileName = "D:/Test/v3d/zzv3d/AF_left_myTest_m.trk";
	vtkSmartPointer<vtkTRKWriter> trkWriter = vtkSmartPointer<vtkTRKWriter>::New();
	trkWriter->SetFileName(trkFileName.c_str());
	trkWriter->SetInputData(reader->GetOutput());

	trkWriter->Write();
	trkWriter->Update();


	return 0;
#endif
}


void ReadNRRD(std::string fileName) {

	vtkNew<vtkTeemNRRDReader> reader;

	reader->SetUseNativeOriginOn();

	reader->SetFileName(fileName.c_str());

	// Check if this is a NRRD file that we can read
	if (!reader->CanReadFile(fileName.c_str()))
	{
		std::cout << ("ReadData: This is not a nrrd file");

		return;
	}

	reader->UpdateInformation();

	// Check type

	if (!(reader->GetPointDataType() == vtkDataSetAttributes::TENSORS))
	{
		std::cout << ("ReadData: MRMLVolumeNode does not match file kind");
		return;
	}

	reader->Update();
	// set volume attributes
	vtkMatrix4x4* mat = reader->GetRasToIjkMatrix();
	// set measurement frame
	vtkMatrix4x4* mat2;
	mat2 = reader->GetMeasurementFrameMatrix();


	double sp[3]{ 0 };
	double orgin[3]{ 0 };
	int dim[3]{ 0 };

	vtkImageData* imgData = reader->GetOutput();
	imgData->GetSpacing(sp);
	imgData->GetOrigin(orgin);
	imgData->GetDimensions(dim);
	std::cout << "spacing:(" << sp[0] << "," << sp[1] << "," << sp[2] << ")" << std::endl;
	std::cout << "orgin:(" << orgin[0] << "," << orgin[1] << "," << orgin[2] << ")" << std::endl;
	std::cout << "dimensions:(" << dim[0] << "," << dim[1] << "," << dim[2] << ")" << std::endl;

}

void ReadNIFTI(std::string fileName,vtkImageData* image) {

	vtkNew<vtkNIFTIImageReader> reader;


	reader->SetFileName(fileName.c_str());

	// Check if this is a NIFTI file that we can read
	if (!reader->CanReadFile(fileName.c_str()))
	{
		std::cout << ("ReadData: This is not a nrrd file");

		return;
	}

	reader->UpdateInformation();


	reader->Update();



	double sp[3]{ 0 };
	double orgin[3]{ 0 };
	int dim[3]{ 0 };

	auto imageData = reader->GetOutput();


	double aFac = reader->GetQFac();
	vtkSmartPointer<vtkMatrix4x4> m_VolumeMatrix = reader->GetQFormMatrix();
	if (m_VolumeMatrix)
	{
		vtkSmartPointer<vtkMatrix4x4> rotateMat = vtkSmartPointer<vtkMatrix4x4>::New();
		rotateMat->SetElement(0, 0, -1);
		rotateMat->SetElement(1, 1, -1);

		vtkMatrix4x4::Multiply4x4(rotateMat, m_VolumeMatrix, m_VolumeMatrix);

		imageData->SetOrigin(
			m_VolumeMatrix->GetElement(0, 3),
			m_VolumeMatrix->GetElement(1, 3),
			m_VolumeMatrix->GetElement(2, 3)
		);

		m_VolumeMatrix->SetElement(0, 3, 0);
		m_VolumeMatrix->SetElement(1, 3, 0);
		m_VolumeMatrix->SetElement(2, 3, 0);

		if (aFac < 0)
		{
			vtkSmartPointer<vtkMatrix4x4> mirrorMat = vtkSmartPointer<vtkMatrix4x4>::New();
			mirrorMat->SetElement(1, 1, -1);
			vtkMatrix4x4::Multiply4x4(mirrorMat, m_VolumeMatrix, m_VolumeMatrix);
		}
	}

	image->DeepCopy(imageData);

}


//----------------------------------------------------------------------------
void ConvertBetweenRASAndLPS(vtkPointSet* inputMesh, vtkPointSet* outputMesh)
{
	vtkNew<vtkTransform> transformRasLps;
	transformRasLps->Scale(-1, -1, 1);
	// vtkTransformPolyDataFilter preserves texture coordinates, while vtkTransformFilter removes them,
	// therefore we must use vtkTransformPolyDataFilter for surface meshes.
	if (inputMesh->IsA("vtkPolyData"))
	{
		vtkNew<vtkTransformPolyDataFilter> transformFilter;
		transformFilter->SetTransform(transformRasLps);
		transformFilter->SetInputData(inputMesh);
		transformFilter->Update();
		outputMesh->ShallowCopy(transformFilter->GetOutput());
	}
	else
	{
		vtkNew<vtkTransformFilter> transformFilter;
		transformFilter->SetTransform(transformRasLps);
		transformFilter->SetInputData(inputMesh);
		transformFilter->TransformAllInputVectorsOn();
		transformFilter->Update();
		outputMesh->ShallowCopy(transformFilter->GetOutput());
	}
}


void LPS2RAS(vtkSmartPointer<vtkImageData> input, vtkMatrix4x4* lpsToIjkMatrix)
{
		unsigned int Dimension = 3;
		vtkMatrix4x4* ijkTolpsMatrix = vtkMatrix4x4::New();

		if (lpsToIjkMatrix == nullptr)
		{
			std::cerr << "ITKWriteVTKImage: rasToIjkMatrix is null" << std::endl;
		}
		else
		{
			vtkMatrix4x4::Invert(lpsToIjkMatrix, ijkTolpsMatrix);
		}

		ijkTolpsMatrix->Transpose();

		double origin[3];
		double mag[3];
		int i;
		for (i = 0; i < 3; i++)
		{
			// normalize vectors
			mag[i] = 0;
			for (int j = 0; j < 3; j++)
			{
				mag[i] += ijkTolpsMatrix->GetElement(i, j) * ijkTolpsMatrix->GetElement(i, j);
			}
			if (mag[i] == 0.0)
			{
				mag[i] = 1;
			}
			mag[i] = sqrt(mag[i]);
		}

		for (i = 0; i < 3; i++)
		{
			int j;
			for (j = 0; j < 3; j++)
			{
				ijkTolpsMatrix->SetElement(i, j, ijkTolpsMatrix->GetElement(i, j) / mag[i]);
			}
		}

		// ITK image direction are in LPS space
		// convert from ijkToRas to ijkToLps
		vtkMatrix4x4* rasToLpsMatrix = vtkMatrix4x4::New();
		rasToLpsMatrix->Identity();
		rasToLpsMatrix->SetElement(0, 0, -1);
		rasToLpsMatrix->SetElement(1, 1, -1);

		vtkMatrix4x4* ijkToLpsMatrix = vtkMatrix4x4::New();
		vtkMatrix4x4::Multiply4x4(ijkTolpsMatrix, rasToLpsMatrix, ijkToLpsMatrix);

		vtkSmartPointer<vtkMatrix3x3> direction = vtkSmartPointer<vtkMatrix3x3>::New();

		for (i = 0; i < Dimension; i++)
		{
			origin[i] = ijkTolpsMatrix->GetElement(3, i);
			int j;
			for (j = 0; j < Dimension; j++)
			{
				if (Dimension == 2)
				{
					//direction[j][i] = (i == j) ? 1. : 0;
					direction->SetElement(i, j, (i == j) ? 1. : 0);
				}
				else
				{
					direction->SetElement(i, j, ijkToLpsMatrix->GetElement(i, j));
					/*direction[j][i] = ijkToLpsMatrix->GetElement(i, j);*/
				}
			}
		}

		rasToLpsMatrix->Delete();
		ijkTolpsMatrix->Delete();
		ijkToLpsMatrix->Delete();

		origin[0] *= -1;
		origin[1] *= -1;

		
		/*input->SetOrigin(origin);
		input->SetSpacing(mag);*/

		vtkNew<vtkImageChangeInformation> ici;
		ici->SetInputData(input);
		ici->SetOutputSpacing(mag);
		ici->SetOutputOrigin(origin);
		ici->Update();
		input->DeepCopy(ici->GetOutput());
		//input->SetDirectionMatrix(direction);
}
