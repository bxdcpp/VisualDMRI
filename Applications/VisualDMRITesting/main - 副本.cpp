#include "vtkSeedTracts.h"

//VTK include
#include <vtkNrrdReader.h>
#include <vtkImageChangeInformation.h>
#include <vtkMatrix4x4.h>
#include <vtkTransform.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkPLYWriter.h>

#include <iostream>
#include <vector>
#include <array>

void CreateTractsForOneSeed(vtkSeedTracts* seed, int thresholdMode,
    double stoppingValue,
    double stoppingCurvature,
    double integrationStepLength,
    double minPathLength,
    double regionSize, double sampleStep,
    int maxNumberOfSeeds,
    int seedSelectedFiducials,double spa[3])
{
	// 设置 world to tensorIJK
	vtkNew<vtkMatrix4x4> mat;
	//volumeNode->GetRASToIJKMatrix(mat.GetPointer());

	vtkNew<vtkMatrix4x4> tensorRASToIJK;
	/*tensorRASToIJK->DeepCopy(mat.GetPointer());*/

	vtkNew<vtkTransform> trans;
	trans->Identity();
	trans->PreMultiply();
	/*trans->SetMatrix(tensorRASToIJK.GetPointer());*/
	// Trans from IJK to RAS
	trans->Inverse();
	// Take into account spacing to compute Scaled IJK
	trans->Scale(1 / spa[0], 1 / spa[1], 1 / spa[2]);
	trans->Inverse();
	seed->SetWorldToTensorScaledIJK(trans.GetPointer());




	//Set tensor Rotation matrix 
	vtkNew<vtkMatrix4x4> tensorRASToIJKRotation;
	tensorRASToIJKRotation->DeepCopy(tensorRASToIJK.GetPointer());

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

	//inputTensorField->GetPointData()->SetScalars(volumeNode->GetImageData()->GetPointData()->GetScalars());





	//6. Extract PolyData in RAS
	vtkSmartPointer<vtkPolyData> outFibers = vtkSmartPointer<vtkPolyData>::New();
	seed->TransformStreamlinesToRASAndAppendToPolyData(outFibers.GetPointer());

	vtkSmartPointer<vtkPLYWriter> write = vtkSmartPointer<vtkPLYWriter>::New();
	std::string fiberFileName = "T:/scene/fiber.ply";
	write->SetFileName(fiberFileName.c_str());
	write->SetInputData(outFibers.Get());
	write->Write();

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
	std::string fileName = "T:/scene/myDWI.nrrd";
    vtkSmartPointer<vtkNrrdReader> nrrdReader = vtkSmartPointer<vtkNrrdReader>::New();
    nrrdReader->SetFileName(fileName.c_str());
	nrrdReader->Update();

	double sp[3]{0};
	vtkImageData* img = nrrdReader->GetOutput();
	if (!img)
	{
		std::cout << "image data is null！" << std::endl;
		return 0;
	}

	img->GetSpacing(sp);
	std::cout << "spacing:(" << sp[0] << "," << sp[1] << "," << sp[2]  << ")";

    vtkNew<vtkImageChangeInformation> ici;
	//1.9531249999999991,1.9531249999999936,3.5100013269750843

    ici->SetOutputSpacing(sp);
    ici->SetInputConnection(nrrdReader->GetOutputPort());
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
		minPathLength, regionSize, sampleStep, maxNumberOfSeeds, seedSelectedFiducials,sp);

	std::cout << "hello visualDMRI" << std::endl;
	return 0;
}