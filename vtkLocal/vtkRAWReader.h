/*=========================================================================

 *	File: vtkRAWReader.h
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
// .NAME vtkRAWReader - read a RAW format file into a VTK structured grid
// .SECTION Description
// Given a RAW format file,  this filter reads the file and covert the 
// structure into vtkImageData. 

#ifndef __vtkRAWReader_h
#define __vtkRAWReader_h

#include "vtkImageAlgorithm.h"
#include "vtkMETAModule.h"
#include "vtkStructuredPoints.h"


class VTK_META_EXPORT vtkRAWReader : public vtkImageAlgorithm
{
public:
  static vtkRAWReader *New();
  vtkTypeMacro(vtkRAWReader,vtkImageAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify file name of the RAW file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  // Directly set the raw file name and dimensionality. 
  void SetRAWFile(char *fileName, int dimX, int dimY, int dimZ);

  // Description:
  // Set the number of dimensions. 
  vtkSetMacro(DimX, int);
  vtkSetMacro(DimY, int);
  vtkSetMacro(DimZ, int);
  
  // Description:
  // Get the number of dimensions.  
  vtkGetMacro(DimX, int);
  vtkGetMacro(DimY, int);
  vtkGetMacro(DimZ, int);
  
  // Description:
  virtual int RequestInformation(
      vtkInformation *, 
      vtkInformationVector **,
      vtkInformationVector *);

protected:
  vtkRAWReader();
  ~vtkRAWReader();
  
  // Name of the input file
  char *FileName;
  // Number of dimensions. 
  int DimX, DimY, DimZ;  
  virtual int RequestData(
      vtkInformation* request,
      vtkInformationVector** inputVector,
      vtkInformationVector* outputVector);
  
private:
  vtkRAWReader(const vtkRAWReader&);    // Not implemented.
  void operator=(const vtkRAWReader&);  // Not implemented.
};

#endif

