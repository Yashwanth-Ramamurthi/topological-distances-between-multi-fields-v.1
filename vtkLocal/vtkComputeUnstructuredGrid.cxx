/*=========================================================================

 *	File: vtkComputeUnstructuredGrid.cxx
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkImageData.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkComputeUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include "vtkObjectFactory.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMath.h"

vtkStandardNewMacro(vtkComputeUnstructuredGrid);

//----------------------------------------------------------------------------
vtkComputeUnstructuredGrid::vtkComputeUnstructuredGrid()
{

}

//----------------------------------------------------------------------------
vtkComputeUnstructuredGrid::~vtkComputeUnstructuredGrid()
{
}

//----------------------------------------------------------------------------
int vtkComputeUnstructuredGrid::FillInputPortInformation(int port, vtkInformation* info)
{
	if(port == 0){
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
	}

  return 1;
}

// vtkComputeUnstructuredGrid - Converts an input dataset consisting of tetrahdera into a vtkUnstructuredGrid
void vtkComputeUnstructuredGrid::ComputeUnstructuredGridFromTetrahedra(vtkFloatArray* pointList, vtkPointData* pointData, vtkIdTypeArray* tetrahedraList){
 	vtkIdType tetrahedraListSize = tetrahedraList->GetNumberOfTuples();
  	vtkIdType pointListSize = pointList->GetNumberOfTuples();
   	output->Allocate(tetrahedraListSize/4);
  	vtkIdType pntIds[4];
	for(int i = 0; i < tetrahedraListSize; i = i + 4){
		pntIds[0] = tetrahedraList->GetValue(i);
                pntIds[1] = tetrahedraList->GetValue(i + 1);
                pntIds[2] = tetrahedraList->GetValue(i + 2);
                pntIds[3] = tetrahedraList->GetValue(i + 3);
		output->InsertNextCell(VTK_TETRA, 4, pntIds);
	}
 	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
 	vtkIdType pnt = 0;
 	points->SetNumberOfPoints(pointListSize/3);
	for(int i = 0; i < pointListSize; i = i + 3){
		points->SetPoint(pnt++, pointList->GetValue(i), pointList->GetValue(i+1), pointList->GetValue(i+2));
	}
	output->SetPoints(points);
 	output->GetPointData()->PassData(pointData);
 	// Avoid keeping extra memory around.
 	output->Squeeze();
}

// vtkComputeUnstructuredGrid - Converts an input dataset consisting of triangles into a vtkUnstructuredGrid
void vtkComputeUnstructuredGrid::ComputeUnstructuredGridFromTriangles(vtkFloatArray* pointList, vtkPointData* pointData, vtkIdTypeArray* triangleList){
 	vtkIdType triangleListSize = triangleList->GetNumberOfTuples();
  	vtkIdType pointListSize = pointList->GetNumberOfTuples();
   	output->Allocate(triangleListSize/3);
  	vtkIdType pntIds[3];
	for(int i = 0; i < triangleListSize; i = i + 3){
		pntIds[0] = triangleList->GetValue(i);
                pntIds[1] = triangleList->GetValue(i + 1);
                pntIds[2] = triangleList->GetValue(i + 2);
		output->InsertNextCell(VTK_TRIANGLE, 3, pntIds);
	}
 	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
 	vtkIdType pnt = 0;
 	points->SetNumberOfPoints(pointListSize/3);
	for(int i = 0; i < pointListSize; i = i + 3){
		points->SetPoint(pnt++, pointList->GetValue(i), pointList->GetValue(i+1), pointList->GetValue(i+2));
	}
	output->SetPoints(points);
 	output->GetPointData()->PassData(pointData);
 	// Avoid keeping extra memory around.
 	output->Squeeze();
}
//----------------------------------------------------------------------------
int vtkComputeUnstructuredGrid::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
	this->output =  vtkUnstructuredGrid::GetData(outputVector);
  return 1;
}

