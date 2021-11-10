/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTRKReader.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkTRKReader.h"


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

vtkStandardNewMacro(vtkTRKReader);

//------------------------------------------------------------------------------
vtkTRKReader::vtkTRKReader()
{
    this->FileName = nullptr;
    this->SetNumberOfInputPorts(0);
}

//------------------------------------------------------------------------------
vtkTRKReader::~vtkTRKReader()
{
    this->SetFileName(nullptr);
}

//------------------------------------------------------------------------------
void vtkTRKReader::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
    os << indent << "FileName: " << (this->FileName ? this->FileName : "(none)") << "\n";
}

//------------------------------------------------------------------------------
int vtkTRKReader::RequestData(
    vtkInformation*, vtkInformationVector**, vtkInformationVector* outputVector)
{
    // Make sure we have a file to read.
    if (!this->FileName)
    {
        vtkErrorMacro("A FileName must be specified.");
        return 0;
    }

  
    // Open the input file.
    m_TrkReader.setFilepath(this->FileName);

    if (!m_TrkReader.open())
    {
        vtkErrorMacro("Error opening file " << this->FileName);
        return 0;
    }
    size_t iTrkTotal = m_TrkReader.getTotalTrkNum();
    // Allocate objects to hold points and vertex cells.
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();

    /// cleanup previous result
    points->Initialize();
    verts->Initialize();

    /// random sampling
    vector<float> cTrk;
    cTrk.reserve(100 * 3);


    vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
    for (int iTrkIdx = 0; iTrkIdx < iTrkTotal; iTrkIdx++)
    {
        cTrk.clear();
        m_TrkReader.readTrack(iTrkIdx, cTrk); ///< read one track from input

        size_t iTotalPnts = cTrk.size() / 3;

        polyLine->Initialize();
        polyLine->GetPointIds()->SetNumberOfIds(iTotalPnts);
        for (int j = 0; j < iTotalPnts; j++)
        {
            vtkIdType iPtIdx = points->InsertNextPoint(cTrk[j * 3], cTrk[j * 3 + 1], cTrk[j * 3 + 2]);
            polyLine->GetPointIds()->SetId(j, iPtIdx);
        }
        verts->InsertNextCell(polyLine);
    }

    // Store the points and cells in the output data object.
    vtkPolyData* pcPolyData = vtkPolyData::GetData(outputVector);

    // Add the points to the dataset
    pcPolyData->SetPoints(points.GetPointer());
    // Add the lines to the dataset
    pcPolyData->SetLines(verts.GetPointer());

    // Add voxel size info (for calcuating the connection table)
    vtkSmartPointer<vtkFloatArray> pcVoxelSize = vtkSmartPointer<vtkFloatArray>::New();
    pcVoxelSize->SetNumberOfComponents(3);
    pcVoxelSize->SetName("voxel_size");
    pcVoxelSize->InsertNextTuple(m_TrkReader.getHeader().voxel_size);
    pcPolyData->GetFieldData()->AddArray(pcVoxelSize);

    return 1;
}
