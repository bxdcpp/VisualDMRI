#include "vtkSeedTracts.h"

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

#include <iostream>
#include <vector>
#include <array>

void ConvertBetweenRASAndLPS(vtkPointSet* inputMesh, vtkPointSet* outputMesh);
void ReadNRRD(std::string fileName);
void seedStreamlineFromPoint(vtkSeedTracts* seed, std::vector<std::array<double, 3>>& seedPoint, double regionSize, double sampleStep);
void CreateTractsForOneSeed(vtkSeedTracts* seed, int thresholdMode,
    double stoppingValue,
    double stoppingCurvature,
    double integrationStepLength,
    double minPathLength,
    double regionSize, double sampleStep,
    int maxNumberOfSeeds,
    int seedSelectedFiducials,double spa[3], vtkMatrix4x4* RASToIJK,vtkImageData* imgData)
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

	std::array<double, 3> points{ -1.160122164827107,-29.925899726775967,60.27479726775956 };
	//std::array<double, 3> points{ 0.39,30.21,59.68 };
	std::vector<std::array<double, 3>> sPoints;
	sPoints.push_back(points);
	seedStreamlineFromPoint(seed, sPoints, regionSize, sampleStep);

	//6. Extract PolyData in RAS
	vtkSmartPointer<vtkPolyData> outFibers = vtkSmartPointer<vtkPolyData>::New();
	seed->TransformStreamlinesToRASAndAppendToPolyData(outFibers.GetPointer());

	

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
	std::string fiberFileName = "W:/scene/fiber.vtk";
	writer->SetFileName(fiberFileName.c_str());
#if VTK_MAJOR_VERSION >= 9
	// version 5.1 is not compatible with earlier Slicer versions and most other software
	writer->SetFileVersion(42);
#endif
	writer->SetFileType( VTK_BINARY);
	
	writer->SetInputData(vtkPointSet::SafeDownCast(outFibers.Get()));
	std::string header = std::string("3D Slicer output. ") + coordinateSytemSpecification;
	writer->SetHeader(header.c_str());
	writer->Write();
	writer->Update();

}

void seedStreamlineFromPoint(vtkSeedTracts* seed, std::vector<std::array<double, 3>>& seedPoint,double regionSize,double sampleStep) 
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
	//0.create vtkSeedTracts
	vtkNew<vtkSeedTracts> seed;
	//1.read nrrd
	std::string fileName = "W:/scene/myDTI.nrrd";
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

	CreateTractsForOneSeed(seed, thresholdMode, stoppingValue, stoppingCurvature, integrationStepLength,
		minPathLength, regionSize, sampleStep, maxNumberOfSeeds, seedSelectedFiducials,sp, mat, ici->GetOutput());

	std::cout << "hello visualDMRI" << std::endl;
	return 0;
}


void ReadNRRD(std::string fileName) {
	
	vtkNew<vtkTeemNRRDReader> reader;

	reader->SetUseNativeOriginOn();

	reader->SetFileName(fileName.c_str());

	// Check if this is a NRRD file that we can read
	if (!reader->CanReadFile(fileName.c_str()))
	{
		std::cout <<("ReadData: This is not a nrrd file");
		
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