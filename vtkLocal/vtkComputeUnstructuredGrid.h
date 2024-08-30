/*=========================================================================

 *	File: vtkComputeUnstructuredGrid.h
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
// .NAME vtkComputeUnstructuredGrid - A variant of vtkSimplicate used to generate an unstructured grid of simplices representing a 3D shape (2-manifold without boundary).
// .SECTION Description
// vtkComputeUnstructuredGrid - Converts an 3D shape (consisting of triangles) into a vtkUnstructuredGrid


#ifndef __vtkComputeUnstructuredGrid_h
#define __vtkComputeUnstructuredGrid_h

#include "vtkMETAModule.h"
#include "vtkUnstructuredGridAlgorithm.h"
#include "vtkIdTypeArray.h"
#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"


class VTK_META_EXPORT vtkComputeUnstructuredGrid : public vtkUnstructuredGridAlgorithm
{
public:
	static vtkComputeUnstructuredGrid *New();
	vtkTypeMacro(vtkComputeUnstructuredGrid,vtkUnstructuredGridAlgorithm);
	
	// vtkComputeUnstructuredGrid - Converts an input dataset consisting of triangles into a vtkUnstructuredGrid
	void ComputeUnstructuredGridFromTriangles(vtkFloatArray* pointList, vtkPointData* pointData, vtkIdTypeArray* triangleList);
	
	// vtkComputeUnstructuredGrid - Converts an input dataset consisting of tetrahdera into a vtkUnstructuredGrid
	void ComputeUnstructuredGridFromTetrahedra(vtkFloatArray* pointList, vtkPointData* pointData, vtkIdTypeArray* tetrahedraList);
	vtkUnstructuredGrid* output;

protected:
	vtkComputeUnstructuredGrid();
	~vtkComputeUnstructuredGrid();

	// Override to specify support for any vtkDataSet input type.
	virtual int FillInputPortInformation(int port, vtkInformation* info);

	virtual int RequestData(
	vtkInformation*,
	vtkInformationVector**,
	vtkInformationVector*);

private:
	  void operator=(const vtkComputeUnstructuredGrid&);  // Not implemented.
};

#endif

