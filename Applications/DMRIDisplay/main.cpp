#include "vtkSeedTracts.h"
#include "vtkTRKReader.h"
#include "vtkTRKWriter.h"
#include "vtkPolyDataColorLinesByOrientation.h"
#include "vtkPolyDataTensorToColor.h"

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
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkLookupTable.h>
#include <vtkProperty.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataReader.h>

#include <iostream>
#include <vector>
#include <array>

//******************************
// 颜色可以设置成三种模式
// Solid Color model
// Mean Qrientation model
// segment Orientation((point_current-point_start)) model
// 
// 显示模式 Line/Tube
// 
// 
//******************************
int main(int argc, char* argv[])
{
#if 0
	int ix = 0;
	std::cout << "hello!" << std::endl;
	vtkSmartPointer<vtkTRKReader> reader1 = vtkSmartPointer<vtkTRKReader>::New();
	std::cout <<reader1->GetReferenceCount()<<"\n";
	//reader1->Delete();
	std::cout << reader1->GetReferenceCount() << "\n";
	std::cout << (reader1.GetPointer() ? 1:0) << "\n";
	std::cout << "modify time:" << reader1->GetMTime() << std::endl;

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	std::cout << actor->GetReferenceCount()<<std::endl;
	actor->Delete();
	std::cout <<"actor:"<< actor->GetReferenceCount() << std::endl;

	vtkSmartPointer<vtkPolyDataWriter> dataWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
	std::cout << dataWriter->GetReferenceCount() << std::endl;

	vtkSmartPointer<vtkPolyDataReader> dataReader = vtkSmartPointer<vtkPolyDataReader>::New();
	std::cout << dataReader->GetReferenceCount() << std::endl;

	vtkSmartPointer<vtkTeemNRRDReader> nrrdreader = vtkSmartPointer<vtkTeemNRRDReader>::New();
	std::cout << nrrdreader->GetReferenceCount() << std::endl;


	return 0;
#else
	///0.read file
	vtkSmartPointer<vtkTRKReader> reader = vtkSmartPointer<vtkTRKReader>::New();
	//const char* fileName = "D:/Test/v3d/zzv3d/AF_left_myTest_m.trk";
	const char* fileName = "C:\\Users\\Bxd\\Documents\\VisualScenes\\20220105152745\\0zcMvK0nStOQDiu8oZcMig.trk";
	reader->SetFileName(fileName);
	reader->Update();


	///
	vtkNew<vtkPolyDataTensorToColor> tensorToColor;
	//是否提取Scalar，根据Scalar 设置颜色
	//tensorToColor->SetExtractScalar(1);
	//tensorToColor->ColorGlyphsByTrace();
	tensorToColor->SetInputConnection(reader->GetOutputPort());

	vtkSmartPointer<vtkPolyDataColorLinesByOrientation>  colorLinesByOrientation = vtkSmartPointer<vtkPolyDataColorLinesByOrientation>::New();
	colorLinesByOrientation->SetColorMode(vtkPolyDataColorLinesByOrientation::colorModeMeanFiberOrientation);
	colorLinesByOrientation->SetInputConnection(tensorToColor->GetOutputPort());

	//管道
	vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
	tubeFilter->SetInputConnection(colorLinesByOrientation->GetOutputPort());
	tubeFilter->SetRadius(0.5);
	tubeFilter->SetNumberOfSides(6);//管道的边数

	vtkNew<vtkLookupTable> colorTable;
	colorTable->SetRampToLinear();
	colorTable->SetTableRange(0, 255);
	colorTable->SetHueRange(0, 1);// 色相
	colorTable->SetSaturationRange(1, 1);//饱和度
	colorTable->SetValueRange(1, 1);
	colorTable->SetAlphaRange(1, 1); // not used
	colorTable->Build();

	

	///1 Create a mapper and actor
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	//mapper->SetInputConnection(reader->GetOutputPort());
	mapper->SetInputConnection(colorLinesByOrientation->GetOutputPort());
	//mapper->SetInputConnection(tubeFilter->GetOutputPort());
	//mapper->SetColorModeToMapScalars();
	mapper->SetScalarRange(0, 255);
	mapper->SetLookupTable(colorTable);

	vtkNew<vtkActor> actor;
	actor->SetMapper(mapper);
	//actor->GetProperty()->SetOpacity(0.2);
	//actor->SetVisibility(0);
	//actor->GetProperty()->SetColor(1, 0, 0);


	///2.Create  a renderer and reder window
	vtkNew<vtkRenderer> renderer;
	vtkNew<vtkRenderWindow> renderWindow;
	renderer->AddActor(actor);
	//renderer->SetBackground()
	renderWindow->AddRenderer(renderer);
	renderWindow->SetWindowName("FiberBundlesWidget");

	/// An interactor
	vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
	renderWindowInteractor->SetRenderWindow(renderWindow);


	/// interactor style.
	vtkNew<vtkInteractorStyleTrackballCamera> style;
	renderWindowInteractor->SetInteractorStyle(style);



	///
	//renderWindow->Render();
	renderWindowInteractor->Initialize();
	renderWindowInteractor->Start();

	return 0;
#endif
}
