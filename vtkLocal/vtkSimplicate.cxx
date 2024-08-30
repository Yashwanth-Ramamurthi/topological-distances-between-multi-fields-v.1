/*=========================================================================

 *	File: vtkSimplicate.cxx
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
#include "vtkIdTypeArray.h"
#include "vtkSimplicate.h"
#include "vtkSmartPointer.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMath.h"
#include "vtkFloatArray.h"


#define CellIndex(i,j,k)   (dim[0] * (dim[1] * (k) + (j)) + (i))
#define CellIndex2(i,j)    (dimA * (j) + (i))

vtkStandardNewMacro(vtkSimplicate);

//----------------------------------------------------------------------------
vtkSimplicate::vtkSimplicate()
{
  this->Scheme = 0;
}

//----------------------------------------------------------------------------
vtkSimplicate::~vtkSimplicate()
{
}

//----------------------------------------------------------------------------
void vtkSimplicate::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  cout << indent << "Simplification scheme: " << this->Scheme << endl;
}

//----------------------------------------------------------------------------
int vtkSimplicate::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
  return 1;
}

//----------------------------------------------------------------------------
int vtkSimplicate::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  // Get input and output data.
  vtkImageData* input = vtkImageData::GetData(inputVector[0]);
  vtkUnstructuredGrid* output = vtkUnstructuredGrid::GetData(outputVector);

  // Skip execution if there is no input geometry.
  vtkIdType numCells = input->GetNumberOfCells();
  vtkIdType numPts = input->GetNumberOfPoints();

  //cout << "numCells - " << numCells << endl;
  //cout << "numPts - " << numPts << endl;

  vtkSmartPointer<vtkIdTypeArray> nbs
      = vtkSmartPointer<vtkIdTypeArray>::New();

  int dim[4];
  vtkIdType pntIds[4];
  
  pntIds[0] = CellIndex(0,0,0);
  
  input->GetDimensions(dim);  
  
  if (numCells < 1 || numPts < 1)
    {
    vtkDebugMacro("No input data.");
    return 1;
    }

  if (input->GetDataDimension() == 3)
    {
    output->Allocate(numCells*6);
    nbs->SetNumberOfComponents(4);
    nbs->SetNumberOfTuples(numCells);
    
    bool render; 
    double p0, p1, p2, p3;
    
    switch (this->Scheme)
      {
      case 0:
        nbs->SetNumberOfTuples(numCells*6);
        for (int k = 0; k < dim[2]-1; k++)
          {
          for (int j = 0; j < dim[1]-1; j++)
            {
            for (vtkIdType i = 0; i < dim[0]-1; i++)
              {
                
                pntIds[0] = CellIndex(i,   j+1, k);
                pntIds[1] = CellIndex(i+1, j+1, k);
                pntIds[2] = CellIndex(i+1, j+1, k+1);
                pntIds[3] = CellIndex(i, j, k);
                
                render = true;
                for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                   {
                   p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                   p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                   p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                   p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                   if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                     {
                        render = false;
                        break;
                     }
                   }
                
               if (render)
                 {
                 output->InsertNextCell(VTK_TETRA, 4, pntIds);
                 }
      
  
               pntIds[0] = CellIndex(i,   j, k);
               pntIds[1] = CellIndex(i+1, j,   k);
               pntIds[2] = CellIndex(i+1, j+1, k+1);
               pntIds[3] = CellIndex(i+1, j+1, k);
                
               render = true;
               for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                  {
                  p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                  p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                  p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                  p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                  if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                    {
                    render = false;
                    break;
                    }
                  }
                
               if (render)
                 {
                 output->InsertNextCell(VTK_TETRA, 4, pntIds);
                 }
           
               pntIds[0] = CellIndex(i,   j,   k);
               pntIds[1] = CellIndex(i+1, j,   k);
               pntIds[2] = CellIndex(i+1, j+1, k+1);
               pntIds[3] = CellIndex(i+1, j,   k+1);
                
               render = true;
               for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                  {
                  p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                  p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                  p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                  p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                  if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                    {
                    render = false;
                    break;
                    }
                  }
                
               if (render)
                 {
                 output->InsertNextCell(VTK_TETRA, 4, pntIds);
                 }
    
               pntIds[0] = CellIndex(i,   j+1, k+1);  // e
               pntIds[1] = CellIndex(i+1, j+1, k+1);  // f
               pntIds[2] = CellIndex(i,   j+1, k);    // a
               pntIds[3] = CellIndex(i,   j,   k);    // c
                
               render = true;
               for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                  {
                  p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                  p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                  p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                  p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                  if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                    {
                    render = false;
                    break;
                    }
                   }
               if (render)
                 {
                 output->InsertNextCell(VTK_TETRA, 4, pntIds);
                 }
                
      
               pntIds[0] = CellIndex(i,   j+1, k+1);  // e
               pntIds[1] = CellIndex(i,   j,   k+1);  // g
               pntIds[2] = CellIndex(i,   j,   k);    // c
               pntIds[3] = CellIndex(i+1, j+1, k+1);  // f
                
               render = true;
               for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                  {
                  p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                  p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                  p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                  p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                  if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                    {
                    render = false;
                    break;
                    }
                  }
                
               if (render)
                 {
                 output->InsertNextCell(VTK_TETRA, 4, pntIds);
                 }
    
               pntIds[0] = CellIndex(i,   j,   k);  // c
               pntIds[1] = CellIndex(i,   j,   k+1);  // g
               pntIds[2] = CellIndex(i+1, j+1, k+1);  // f
               pntIds[3] = CellIndex(i+1, j,   k+1);  // h
                
               render = true;
               for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                  {
                  p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                  p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                  p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                  p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                  if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                    {
                    render = false;
                    break;
                    }
                  }
                if (render)
                  {
                  output->InsertNextCell(VTK_TETRA, 4, pntIds);
                  }
              }
            }
          }
        break;
        
      case 1:
        nbs->SetNumberOfTuples(numCells*5);
        for (int k = 0; k < dim[2]-1; k++)
          {
          for (int j = 0; j < dim[1]-1; j++)
            {
            for (vtkIdType i = 0; i < dim[0]-1; i++)
              {
              if ((i+j+k) % 2)
                {
                pntIds[0] = CellIndex(i,   j,   k);
                pntIds[1] = CellIndex(i,   j+1, k);
                pntIds[2] = CellIndex(i+1, j+1, k);
                pntIds[3] = CellIndex(i,   j+1, k+1);
                render = true;
                for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                   {
                   p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                   p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                   p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                   p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                   if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                     {
                     render = false;
                     break;
                     }
                   }
                if (render)
                  {
                  output->InsertNextCell(VTK_TETRA, 4, pntIds);
                  }
      
                pntIds[0] = CellIndex(i,   j+1, k+1);
                pntIds[1] = CellIndex(i+1, j+1, k+1);
                pntIds[2] = CellIndex(i+1, j+1, k  );
                pntIds[3] = CellIndex(i+1, j,   k+1);
                render = true;
                for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                   {
                   p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                   p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                   p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                   p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                   if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                     {
                     render = false;
                     break;
                     }
                   }
                if (render)
                  {
                  output->InsertNextCell(VTK_TETRA, 4, pntIds);
                  }
                  
                pntIds[0] = CellIndex(i,   j,   k  );
                pntIds[1] = CellIndex(i+1, j,   k  );
                pntIds[2] = CellIndex(i+1, j+1, k  );
                pntIds[3] = CellIndex(i+1, j,   k+1);
               
                render = true;
                for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                   {
                   p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                   p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                   p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                   p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                   if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                     {
                     render = false;
                     break;
                     }
                   }
                 if (render)
                   {
                   output->InsertNextCell(VTK_TETRA, 4, pntIds);
                   }

                 pntIds[0] = CellIndex(i,   j,   k);
                 pntIds[1] = CellIndex(i,   j,   k+1);
                 pntIds[2] = CellIndex(i+1, j,   k+1);
                 pntIds[3] = CellIndex(i,   j+1, k+1);
                    
                render = true;
                for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                   {
                   p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                   p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                   p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                   p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                   if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                     {
                     render = false;
                     break;
                     }
                   }
                 if (render)
                   {
                   output->InsertNextCell(VTK_TETRA, 4, pntIds);
                   }

                 pntIds[0] = CellIndex(i,   j,   k  );
                 pntIds[1] = CellIndex(i+1, j,   k+1);
                 pntIds[2] = CellIndex(i,   j+1, k+1);
                 pntIds[3] = CellIndex(i+1, j+1, k  );
                render = true;
                for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                   {
                   p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                   p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                   p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                   p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                   if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                     {
                     render = false;
                     break;
                     }
                   }
                 if (render)
                   {
                   output->InsertNextCell(VTK_TETRA, 4, pntIds);
                   }    
                }
              else
                {
                pntIds[0] = CellIndex(i,   j,   k);
                pntIds[1] = CellIndex(i,   j+1, k);
                pntIds[2] = CellIndex(i+1, j,   k);
                pntIds[3] = CellIndex(i,   j,   k+1);
                
               render = true;
               for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                  {
                  p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                  p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                  p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                  p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                  if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                    {
                    render = false;
                    break;
                    }
                   }
                
                if (render)
                  {
                  output->InsertNextCell(VTK_TETRA, 4, pntIds);
                  }
      
                pntIds[0] = CellIndex(i,   j,   k+1);
                pntIds[1] = CellIndex(i+1, j+1, k+1);
                pntIds[2] = CellIndex(i+1, j,   k  );
                pntIds[3] = CellIndex(i+1, j,   k+1);
                render = true;
                for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                   {
                   p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                   p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                   p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                   p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                   if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                     {
                     render = false;
                     break;
                     }
                   }
                if (render)
                  {
                  output->InsertNextCell(VTK_TETRA, 4, pntIds);
                  }
    
                pntIds[0] = CellIndex(i+1, j+1, k+1);
                pntIds[1] = CellIndex(i+1, j,   k  );
                pntIds[2] = CellIndex(i+1, j+1, k  );
                pntIds[3] = CellIndex(i,   j+1, k  );
                render = true;
                for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                   {
                   p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                   p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                   p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                   p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                   if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                     {
                     render = false;
                     break;
                     }
                   }
                if (render)
                  {
                  output->InsertNextCell(VTK_TETRA, 4, pntIds);
                  }
    
                pntIds[0] = CellIndex(i,   j+1, k);
                pntIds[1] = CellIndex(i+1, j+1, k+1);
                pntIds[2] = CellIndex(i,   j+1, k+1);
                pntIds[3] = CellIndex(i,   j,   k+1);
                render = true;
                for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                   {
                   p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                   p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                   p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                   p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                   if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                     {
                     render = false;
                     break;
                     }
                   }
                if (render)
                  {
                  output->InsertNextCell(VTK_TETRA, 4, pntIds);
                  }
      
                pntIds[0] = CellIndex(i,   j,   k+1);
                pntIds[1] = CellIndex(i+1, j,   k  );
                pntIds[2] = CellIndex(i+1, j+1, k+1);
                pntIds[3] = CellIndex(i,   j+1, k  );
                render = true;
                for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
                   {
                   p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
                   p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
                   p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
                   p3 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[3],0);
                   if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2)
                       || vtkMath::IsNan(p3))
                     {
                     render = false;
                     break;
                     }
                   }
                if (render)
                  {
                  output->InsertNextCell(VTK_TETRA, 4, pntIds);
                  }
                }
              }
            }
          }
        break;

      default:
        vtkWarningMacro("Unsupported 3D schame, no output.");
      }
    }
  else
    {
    int dimA, dimB;
    if (dim[0] == 1)
      {
      dimA = dim[1];
      dimB = dim[2];
      }
    if (dim[1] == 1)
      {
      dimA = dim[0];
      dimB = dim[2];
      }
    if (dim[2] == 1)
      {
      dimA = dim[0];
      dimB = dim[1];
      }
    output->Allocate(numCells*2);
    nbs->SetNumberOfComponents(3);
    nbs->SetNumberOfTuples(numCells*2);

    bool render;
    double p0,p1,p2;
    for (int j = 0; j < dimB-1; j++)
      {
      for (vtkIdType i = 0; i < dimA-1; i++)
        {
        if ((i+j) % 2)
          {
          pntIds[0] = CellIndex2(i,   j+1);
          pntIds[1] = CellIndex2(i,   j  );
          pntIds[2] = CellIndex2(i+1, j  );
          
          render = true;
          for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
             {
             p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
             p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
             p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
             if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2))
               {
               render = false;
               break;
                }
             }
          if (render)
            {
            output->InsertNextCell(VTK_TRIANGLE, 3, pntIds);
            }
          
          pntIds[0] = CellIndex2(i,   j+1);
          pntIds[1] = CellIndex2(i+1, j  );
          pntIds[2] = CellIndex2(i+1, j+1);
          
          render = true;
          for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
             {
             p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
             p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
             p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
             if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2))
               {
               render = false;
               break;
                }
             }
          if (render)
            {
            output->InsertNextCell(VTK_TRIANGLE, 3, pntIds);
            }
          }
        else
          {
          pntIds[0] = CellIndex2(i,   j  );
          pntIds[1] = CellIndex2(i+1, j  );
          pntIds[2] = CellIndex2(i+1, j+1);
          render = true;
          for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
             {
             p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
             p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
             p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
             if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2))
               {
               render = false;
               break;
                }
             }
          if (render)
            {
            output->InsertNextCell(VTK_TRIANGLE, 3, pntIds);
            }

          pntIds[0] = CellIndex2(i,   j  );
          pntIds[1] = CellIndex2(i,   j+1);
          pntIds[2] = CellIndex2(i+1, j+1);
          render = true;
          for (int m = 0; m < input->GetPointData()->GetNumberOfArrays();m++)
             {
             p0 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[0],0);
             p1 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[1],0);
             p2 = input->GetPointData()->GetArray(m)->GetComponent(pntIds[2],0);
             if (vtkMath::IsNan(p0) || vtkMath::IsNan(p1) || vtkMath::IsNan(p2))
               {
               render = false;
               break;
                }
             }
          if (render)
            {
            output->InsertNextCell(VTK_TRIANGLE, 3, pntIds);
            }
          }
        }
      }
    }
  
  // Explicitly construct points.
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkIdType pnt = 0;
  
  //std::cout << "dim0: " << dim[0] << "  dim1: " << dim[1] << "  dim2: " << dim[2] << std::endl;
  points->SetNumberOfPoints(dim[0]*dim[1]*dim[2]);

  for (int k = 0; k < dim[2]; k++)
    {
    for (int j = 0; j < dim[1]; j++)
      {
      for (int i = 0; i < dim[0]; i++)
        {
        points->SetPoint(pnt++, i, j, k);
        }
      }
    }
  // Store the new set of points in the output.
  output->SetPoints(points);

  // Pass point data through because we still have the same number
  // of points.
  output->GetPointData()->PassData(input->GetPointData());
  // Avoid keeping extra memory around.
  output->Squeeze();

  return 1;
}

