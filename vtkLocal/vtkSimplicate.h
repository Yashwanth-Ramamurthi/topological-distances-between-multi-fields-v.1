/*=========================================================================

 *	File: vtkSimplicate.h
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
// .NAME vtkSimplicate - Generate an unstructured grid of simplices.
// .SECTION Description
// vtkSimplicate converts an input dataset into an unstructured grid
// of simplices.  The filter operates in 2D or 3D mode.  For 3D
// voxels, a flag determines whether voxels are subdivided into 
// five or six tets. If the scalar value of the current cell is defined
// Nan, then this cell will be discarded. 


#ifndef __vtkSimplicate_h
#define __vtkSimplicate_h

#include "vtkMETAModule.h"
#include "vtkUnstructuredGridAlgorithm.h"
//#include <CoinPackedVector.hpp>


class VTK_META_EXPORT vtkSimplicate : public vtkUnstructuredGridAlgorithm
{
public:
  static vtkSimplicate *New();
  vtkTypeMacro(vtkSimplicate,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Select subdivision scheme for 3D input.
  // 0 = 6-tet subdivision, 1 = 5-tet subdivision.
  vtkSetMacro(Scheme, int);
  vtkGetMacro(Scheme, int);

protected:
  vtkSimplicate();
  ~vtkSimplicate();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);

  // Main implementation.
  virtual int RequestData(
      vtkInformation*,
      vtkInformationVector**,
      vtkInformationVector*);

  // Subdivision scheme
  int Scheme;

private:
  vtkSimplicate(const vtkSimplicate&);  // Not implemented.
  void operator=(const vtkSimplicate&);  // Not implemented.
};

#endif

