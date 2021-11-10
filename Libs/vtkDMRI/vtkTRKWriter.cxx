
#include "vtkTRKWriter.h"


//VTK include 
#include <vtkCellArray.h>
#include <vtkObjectFactory.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtksys/FStream.hxx>
#include <vtkObjectFactory.h>
#include <vtkPolyLine.h>
#include <vtkFloatArray.h>
#include <vtkFieldData.h>
#include <vtkInformation.h>

vtkStandardNewMacro(vtkTRKWriter);

//------------------------------------------------------------------------------
vtkTRKWriter::vtkTRKWriter()
{
	this->FileName = nullptr;
	this->VoxelOrder[0] = 'R';
	this->VoxelOrder[1] = 'A';
	this->VoxelOrder[2] = 'S';
	this->VoxelOrder[3] = '\0';

	//set trk header default value
	m_TrkHeader.voxel_size[0] = 2.50;
	m_TrkHeader.voxel_size[1] = 2.50;
	m_TrkHeader.voxel_size[2] = 2.50;


	/*this->SetDebug(true);
	this->SetGlobalWarningDisplay(true);*/
}

//------------------------------------------------------------------------------
vtkTRKWriter::~vtkTRKWriter()
{
	delete[] this->FileName;
}

//------------------------------------------------------------------------------
void vtkTRKWriter::WriteData()
{
	vtkIdType i,  idx;
	vtkPoints* inPts;
	vtkCellArray* verts;

	vtkPolyData* input = this->GetInput();

	//Get input and check data
	verts = input->GetLines(); 
	inPts = input->GetPoints();

	if (inPts == nullptr || verts == nullptr)
	{
		vtkErrorMacro(<< "No data to write!");
		return;
	}

	vtkIdType numLines = verts->GetNumberOfCells();
	vtkIdType numPts = inPts->GetNumberOfPoints();

	vtkDebugMacro(<< numLines);
	vtkDebugMacro(<< numPts);

	std::string wFileName;
	wFileName.assign(this->FileName, this->FileName + strlen(this->FileName));
	/// create writer and create an empty new file
	TrkFileWriter cWriter(wFileName);
	strcpy_s(m_TrkHeader.voxel_order, VoxelOrder);
	cWriter.copyHeader(m_TrkHeader);
	if (!cWriter.create())
		return;

	double pointData[3]{ 0 };
	std::vector<float> cTrk;
	for (idx = 0;idx < numLines;idx++)
	{
		vtkNew<vtkIdList> ptsId;
		ptsId->Initialize();
		verts->GetNextCell(ptsId);
		//verts->GetCell(idx, ptsId);//get point id
		vtkIdType ptsNumberofLine = ptsId->GetNumberOfIds();//get point num

		cTrk.clear();
		for (i = 0; i < ptsNumberofLine; i++)
		{
			//vtkDebugMacro(<< ptsId->GetId(i));
			inPts->GetPoint(ptsId->GetId(i), pointData);
			cTrk.push_back(pointData[0]);
			cTrk.push_back(pointData[1]);
			cTrk.push_back(pointData[2]);
	
		}
		cWriter.appendTrack(cTrk);
	}

	cWriter.close();

	//ptsId->Delete();
	//input->GetPolys()

}

//------------------------------------------------------------------------------
int vtkTRKWriter::FillInputPortInformation(int port, vtkInformation* info)
{
	info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
	return 1;
}

//------------------------------------------------------------------------------
void vtkTRKWriter::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Storing order of the original image data:";
	os << this->VoxelOrder << "\n";
}

//------------------------------------------------------------------------------
vtkPolyData* vtkTRKWriter::GetInput()
{
	return vtkPolyData::SafeDownCast(this->Superclass::GetInput());
}

//------------------------------------------------------------------------------
vtkPolyData* vtkTRKWriter::GetInput(int port)
{
	return vtkPolyData::SafeDownCast(this->Superclass::GetInput(port));
}


