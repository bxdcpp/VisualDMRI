/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkTRKReader.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

	 This software is distributed WITHOUT ANY WARRANTY; without even
	 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
	 PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkTRKReader
 * @brief   read Stanford University TRA polygonal file format
 *
* vtkPLYReader is a source object that reads polygonal data in
 * Stanford TRK file format (see
 * http://www.trackvis.org/docs/?subsect=fileformat).
 *
 */

#ifndef vtkTRKReader_h
#define vtkTRKReader_h

#include <vtkIOLegacyModule.h> // For export macro
#include <vtkPolyDataAlgorithm.h>


#include "vtkDMRIConfigure.h"
#include "defs.h"
#include "trkfileio.h"

class vtkDMRI_EXPORT vtkTRKReader : public vtkPolyDataAlgorithm
{
public:
	static vtkTRKReader* New();
	vtkTypeMacro(vtkTRKReader, vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent) override;

	//@{
	/**
	 * Set/Get the name of the file from which to read points.
	 */
	vtkSetStringMacro(FileName);
	vtkGetStringMacro(FileName);
	//@}

	ADD_CLASS_FIELD_NOSETTER(TrkFileReader, TrkReader, getTrkReader);
	
protected:
	vtkTRKReader();
	~vtkTRKReader() override;

	char* FileName;

	int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
	vtkTRKReader(const vtkTRKReader&) = delete;
	void operator=(const vtkTRKReader&) = delete;
};

#endif
