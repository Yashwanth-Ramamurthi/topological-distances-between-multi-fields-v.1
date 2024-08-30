/*=========================================================================

 *	File: vtkMultiDimensinoalReebGraph.cxx
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
#include "vtkAdjacentVertexIterator.h"
#include "vtkBitArray.h"
#include "vtkBoostConnectedComponents.h"
#include "vtkDataSetAttributes.h"
#include "vtkFloatArray.h"
#include "vtkGraphEdge.h"
#include "vtkIdTypeArray.h"
#include "vtkIntArray.h"
#include "vtkInformation.h"
#include "vtkMultiDimensionalReebGraph.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkMutableUndirectedGraph.h"
#include "vtkMutableGraphHelper.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkTreeBFSIterator.h"
#include "vtkVoidArray.h"
#include <queue>
#include <stack>

vtkStandardNewMacro(vtkMultiDimensionalReebGraph);

#define VTK_CREATE(type,name) \
vtkSmartPointer<type> name = vtkSmartPointer<type>::New()


//-----------------------------------------------------------------------
vtkMultiDimensionalReebGraph::vtkMultiDimensionalReebGraph()
{
  this->SetNumberOfOutputPorts(1);
  this->Fields = vtkstd::vector<char*>();
  this->MDRG=NULL;
}

vtkMultiDimensionalReebGraph::~vtkMultiDimensionalReebGraph()
{
  while (this->Fields.size() > 0)
       {
       delete [] this->Fields.back();
       this->Fields.pop_back();
       }
  this->Fields.clear();
}

void vtkMultiDimensionalReebGraph::AddField(char *nm)
{
  char *f = new char[strlen(nm) + 1];
  strcpy(f, nm);
  this->Fields.push_back(f);
}

void vtkMultiDimensionalReebGraph::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

int vtkMultiDimensionalReebGraph::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkGraph");
    return 1;
    }
  return 0;
}

//----------------------------------------------------------------------------------------
//Extract MDRG
vtkTree *vtkMultiDimensionalReebGraph::GetMultiDimensionalReebGraph()
{ 
  if (!this->MDRG)
    {
    vtkErrorMacro("Input has not been processed, no graph available: First Extract the MDRG..");
    return 0;
    }
  return MDRG;
}

//Extract node (which is a ReebGraph) of the MDRG
vtkSmartPointer<vtkExtractReebGraph> vtkMultiDimensionalReebGraph::GetTreeNode(vtkIdType i)
{
	//cout<<"\nNumber of graphs : "<<ReebGraphList.size()<<"\n\n";
  return this->ReebGraphList[i];
}

//Extract Jacobi and Boundary nodes
vtkSmartPointer<vtkIdTypeArray> vtkMultiDimensionalReebGraph::ExtractJacobiNodes()
{
  VTK_CREATE(vtkIdTypeArray, nodeIds);
  for (vtkIdType i = 0; i < this->MDRG->GetNumberOfVertices(); i++)
     {
     if (this->MDRG->IsLeaf(i))
       {
       vtkSmartPointer<vtkExtractReebGraph> ect=GetTreeNode(i);
       vtkSmartPointer<vtkIdTypeArray> arr=ect->GetCriticalNodes();
  
       for (vtkIdType i=0; i<arr->GetNumberOfTuples();  i++)
          {
          nodeIds->InsertNextValue(arr->GetValue(i));
          }
       }
      }
  return nodeIds;
}
//-----------------------------------------------------------------------------------------
void vtkMultiDimensionalReebGraph::CreateMultiDimensionalReebGraph(vtkGraph *inputGraph)
{
  if (!inputGraph)
    {
    vtkErrorMacro("A vtkGraph object is required.");
    } 
  VTK_CREATE(vtkMutableDirectedGraph, dirGraph);	
  
  //Data-structures needed for creating MDRG
  std::queue<vtkGraph*> queueSG;
  std::queue<vtkSmartPointer<vtkExtractReebGraph> >queueRG;
  queueSG.push(inputGraph);

  std::queue<vtkIdType> queueIds;
  vtkIdType id;
  id = dirGraph->AddVertex();
  queueIds.push(id);
 //cout<<"\nGraph id : "<<id;
  //Creating MultiDimensional Reeb Graph
  //This is BFS without recursion
  for (vtkIdType f=0; f<Fields.size(); f++)
     {
     while (!queueSG.empty())
          {
          vtkGraph *Graph=queueSG.front();
          queueSG.pop();
          VTK_CREATE(vtkExtractReebGraph, erg0);
          vtkIdType success= erg0->ExtractReebGraph(Graph, this->Fields[f]);
          queueRG.push(erg0);
          this->ReebGraphList.push_back(erg0);
          }
     while (!queueRG.empty() &&  f<Fields.size()-1)
          {
          vtkSmartPointer<vtkExtractReebGraph> erg1=queueRG.front();
          queueRG.pop();
          vtkReebGraph *RG=erg1->GetReebGraph();

          vtkIdType S = queueIds.front();
          queueIds.pop();

         for (vtkIdType i = 0; i < RG->GetNumberOfVertices(); i++)
            {
            vtkGraph *SG=erg1->GetSubGraph(i);
            queueSG.push(SG);
            vtkIdType T = dirGraph->AddVertex();
 //cout<<"\nGraph id : "<<T<< "Field : "<<f<<"\n\n";
            queueIds.push(T);
            dirGraph->AddEdge(S, T);
            }
          }
      }
  VTK_CREATE(vtkTree, outputTree);  
  bool success = outputTree->CheckedShallowCopy(dirGraph);
  this->MDRG=outputTree;
}

//----------------------------------------------------------------------------
int vtkMultiDimensionalReebGraph::RequestData(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, 
  vtkInformationVector* outputVector)
{
  // Ensure we have valid inputs ...
  vtkGraph* const inputGraph = vtkGraph::GetData(inputVector[0]);
  vtkGraph* const outputGraph = vtkGraph::GetData(outputVector);
  
  if (!inputGraph)
    {
    vtkErrorMacro("A vtkGraph object is required.");
    return 0;
    }

  vtkIdType nrVertices = inputGraph->GetNumberOfVertices();
  if (!nrVertices)
    {
    vtkErrorMacro("Empty graph, no vertices.");
    return 0;
    }

  if (!this->Fields[0])
    {
    vtkErrorMacro("No field specified.");
    return 0;
    }

  //Setting JCN-ids
  VTK_CREATE(vtkIdTypeArray, nodeIds);
  nodeIds->SetName("jcnids");
  nodeIds->SetNumberOfValues(nrVertices);
  for (vtkIdType i = 0; i < nrVertices; i++)
     {
     nodeIds->SetValue(i, i);//copying the origical JCN-ids!	
     double point[3];
     inputGraph->GetPoint(i,point);    
     }
  inputGraph->GetVertexData()->AddArray(nodeIds);//Adding JCN-ids in Input Graph
  CreateMultiDimensionalReebGraph(inputGraph);
  vtkSmartPointer<vtkIdTypeArray> jacobiNodes=ExtractJacobiNodes();

  //Need to set Flag 1: for Jacobi nodes
  VTK_CREATE(vtkIntArray, jacobiFlag);
  jacobiFlag->SetNumberOfTuples(nrVertices);
  jacobiFlag->SetName("jacobi-flag");

  //Initializing Jacobi Set Flags to 0
  for (int i = 0; i < nrVertices; i++)
     {
     jacobiFlag->SetValue(i,0);
     }
  for (vtkIdType i=0; i<jacobiNodes->GetNumberOfTuples(); i++)
     {
     jacobiFlag->SetValue(jacobiNodes->GetComponent(i,0),1);
     }

  //Memory allocate
  vtkSmartPointer<vtkBitArray> value = vtkSmartPointer<vtkBitArray>::New();
  value->SetNumberOfValues(nrVertices);
  value->SetName("MDRG");

  //Setting field-values for Jacobi-nodes
  for (int i = 0; i < nrVertices; i++)
     {
     if (jacobiFlag->GetValue(i)==1)
       {
       value->SetValue(i,1);//jcacobi
       }
     else
       {
       value->SetValue(i,0);//none jacobi		
       }
     }
 
  inputGraph->GetVertexData()->RemoveArray("jcnids");
  vtkDataSetAttributes *inPD   = inputGraph->GetVertexData(), 
  *outPD  = outputGraph->GetVertexData();
  outputGraph->ShallowCopy(inputGraph);
  outPD->CopyAllocate(inPD, nrVertices, 10);	
  for (vtkIdType i = 0; i < nrVertices; i++)
     {
     outPD->CopyData(inputGraph->GetVertexData(), i, i);
     }
      
  //Added for coloring Jacobi Nodes
  outPD->AddArray(value);   
  outputGraph->Squeeze();  
  return 1;
}
//--------------------------------------------------------------------------------------

