/**
 * @class   vtkTRKWriter
 * @brief   write Stanford TRK file format
 *
 * vtkTRKWriter writes polygonal data in Stanford University PLY format
 * (see http://www.trackvis.org/docs/?subsect=fileformat).
 * @warning
 *
 * @sa
 * vtkTRKReader
 */


#ifndef vtkTRKWriter_h
#define vtkTRKWriter_h


#include "vtkDMRIConfigure.h"
#include "trkfileio.h"

#include <vtkIOLegacyModule.h> // For export macro
#include <vtkPolyDataAlgorithm.h>
#include <vtkWriter.h>


#include <string>

class vtkDMRI_EXPORT vtkTRKWriter : public vtkWriter
{
public:
	static vtkTRKWriter* New();
	vtkTypeMacro(vtkTRKWriter, vtkWriter);
	void PrintSelf(ostream& os, vtkIndent indent) override;

	//@{
 /**
  * Get the input to this writer.
  */
	vtkPolyData* GetInput();
	vtkPolyData* GetInput(int port);
	//@}

	//@{
	/**
	 * Specify file name of vtk polygon data file to write.
	 */
	vtkSetStringMacro(FileName);
	vtkGetStringMacro(FileName);
	//@}

	 //@{
  /**
   * Storing order of the original image data. Explained here.
   */
	vtkSetVector3Macro(VoxelOrder, char);
	vtkGetVector3Macro(VoxelOrder, char);
	//@}

	
protected:
	vtkTRKWriter();
	~vtkTRKWriter() override;
	void WriteData() override;
	int FillInputPortInformation(int port, vtkInformation* info) override;


	char* FileName;
	char VoxelOrder[4];

	ADD_CLASS_FIELD_NOSETTER(TrkFileHeader, TrkHeader, GetHeader)     ///< track file header        

private:
	vtkTRKWriter(const vtkTRKWriter&) = delete;
	void operator=(const vtkTRKWriter&) = delete;
};

#endif
