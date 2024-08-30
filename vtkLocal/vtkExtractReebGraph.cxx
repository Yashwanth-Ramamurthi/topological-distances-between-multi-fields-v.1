/*=========================================================================

 *	File: vtkExtractReebGraph.cxx
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
#include "vtkAdjacentVertexIterator.h"
#include "vtkBoostConnectedComponents.h"
#include "vtkBitArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkExtractReebGraph.h"
#include "vtkFloatArray.h"
#include "vtkGraphEdge.h"
#include "vtkIdTypeArray.h"
#include "vtkInformation.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkMutableUndirectedGraph.h"
#include "vtkMutableGraphHelper.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include <queue>
#include <stack>
#include <boost/unordered_map.hpp>

vtkStandardNewMacro(vtkExtractReebGraph);

#define VTK_CREATE(type,name) \
   vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

//------------------------------------------------------------
vtkExtractReebGraph::vtkExtractReebGraph()
{
  this->FieldName = NULL;
  this->UF = vtkSmartPointer<vtkIdTypeArray>::New();
  this->Class = vtkSmartPointer<vtkIdTypeArray>::New();
  this->JCNIndex = vtkSmartPointer<vtkIdTypeArray>::New();
  this->Start = vtkSmartPointer<vtkIdTypeArray>::New();
  this->Size = vtkSmartPointer<vtkIdTypeArray>::New();
  this->Graph = NULL;
  this->ReebGraph = NULL;
}

vtkExtractReebGraph::~vtkExtractReebGraph()
{
  if (this->FieldName)
    {
    delete [] this->FieldName;
    }
}

void vtkExtractReebGraph::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

int vtkExtractReebGraph::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0)
    {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkGraph");
    return 1;
    }
    
  return 0;
}
//-----------------------------------------------------------------------------
// Find the root of a value stored in the union-find
// structure, performing path compression as we go.
vtkIdType vtkExtractReebGraph::Root(vtkIdType x)
{
  vtkIdType entry = this->UF->GetValue(x);
  
  if (entry < 0) //critical node
    {
    return x;
    }
  else
    {
    vtkIdType root = this->Root(entry);
    this->UF->SetValue(x, root); //needed for Path Compression
    return root;
    }
}

// Find if x and y are in the same component
bool  vtkExtractReebGraph::Find(vtkIdType x, vtkIdType y)
{
  vtkIdType root_x = this->Root(x);
  vtkIdType root_y = this->Root(y);

  return (root_x==root_y);
}

//Merge components containing i and nb, if they are not in a same component
void vtkExtractReebGraph::Union(vtkIdType i, vtkIdType nb)
{
  vtkIdType root_i  = this->Root(i);
  vtkIdType root_nb = this->Root(nb);
  if (root_i != root_nb)//not already in the same component
    {
    this->UF->SetValue(root_i, root_nb);//make in same component
    }
}

void vtkExtractReebGraph::CreateUnionFindStructure(vtkGraph *inputGraph, char *FieldName)
{
  vtkDataSetAttributes *inPD   = inputGraph->GetVertexData();
  vtkDataArray *property = inPD->GetArray(FieldName); 
  vtkIdType nrVertices = inputGraph->GetNumberOfVertices();
  this->UF->SetNumberOfTuples(nrVertices);
  VTK_CREATE(vtkBitArray, seen);
  seen->SetNumberOfTuples(nrVertices);

  for (vtkIdType i = 0; i < nrVertices; i++)
     {
     this->UF->SetValue(i, -1);
     seen->SetValue(i, 0);
     }

  double value;
  vtkIdType nb, root_nb, root_i;

  VTK_CREATE(vtkAdjacentVertexIterator, adj);
  // Look at each vertex, and merge adjacent vertices where the specified
  // field components have the same value.
  for (vtkIdType i = 0; i < nrVertices; i++)
     {
     value = property->GetComponent(i, 0);
     inputGraph->GetAdjacentVertices(i, adj);

     while (adj->HasNext())//for nodes in a connected component
          {
          nb = adj->Next();
          // If neighbours have same value for this property, they
          // belong to the cluster of nodes in the contour tree:
          // merge classes.
          
          //If having same value
          if (!seen->GetValue(nb) && property->GetComponent(nb, 0) == value)
            {
            Union(i,nb);
            }
          } // for each neighbour
     seen->SetValue(i, 1);
     } // for each vertex
}

//---------------------------------------------------------------------------------------------------------------
vtkIdType vtkExtractReebGraph::GetSizeOfSubGraph(vtkIdType groupId)
{
  if (!this->Graph)
    {
    vtkErrorMacro("Input has not been processed, no graph available: First Extract the Reeb Graph..");
    return 0;
    }

  vtkIdType size = this->Size->GetValue(groupId);
  return size;
}

vtkSmartPointer<vtkIdTypeArray> vtkExtractReebGraph::GetSubGraphJCNIds(vtkIdType groupId)
{
  if (!this->Graph)
    {
    vtkErrorMacro("Input has not been processed, no graph available: First Extract the Reeb Graph..");
    return 0;
    }

  vtkIdType base = this->Start->GetValue(groupId);
  vtkIdType size = this->Size->GetValue(groupId);
  VTK_CREATE(vtkIdTypeArray, nodeIds);
  nodeIds->SetName("jcnids");
  nodeIds->SetNumberOfValues(size);
  for (vtkIdType i = 0; i < size ; i++)
     {
     vtkIdType v = this->JCNIndex->GetValue(base + i);
     nodeIds->SetValue(i, v);
     }
  return nodeIds;
}


//----------------------------------------------------------------------------
// This function returns a Contour (connected component of a Levet set) 
// corresponding to a node in JCN
vtkGraph *vtkExtractReebGraph::GetSubGraph(vtkIdType groupId)
{
  VTK_CREATE(vtkMutableGraphHelper, builder);
  VTK_CREATE(vtkMutableUndirectedGraph, undir);
  builder->SetGraph(undir);
  
  VTK_CREATE(vtkIdTypeArray, dest);
  vtkDataSetAttributes *data = undir->GetVertexData();
  //Setting the pointer of vertex-data of 'undir' in data

  VTK_CREATE(vtkMutableUndirectedGraph, output);
  vtkIdType base, u, v;

  if (!this->Graph)
    {
    vtkErrorMacro("Input has not been processed, no graph available: First Extract the Reeb Graph..");
    return 0;
    }

  vtkIdType nrVertices = this->Graph->GetNumberOfVertices();

  data->CopyAllocate(this->Graph->GetVertexData(), 
                     this->Size->GetValue(groupId), 
                     10);

  //'dest': an array of same size as parent graph and map node-ids of sub-graph 
  dest->SetNumberOfValues(nrVertices);
  for (vtkIdType i = 0; i < nrVertices; i++)
    {
    dest->SetValue(i, -1);
    }

  //Copying jcnids of the input-graph
  VTK_CREATE(vtkIdTypeArray, nodeIds);
  nodeIds=GetSubGraphJCNIds(groupId);//Need to copy the original JCN-ids, here!!

  //Set data to subgraph.
  base = this->Start->GetValue(groupId);
  vtkIdType vjcnid;
  for (vtkIdType i = 0; i < this->Size->GetValue(groupId); i++)
     {
     builder->AddVertex();
     //get Parent Graph node-ids
     v = this->Class->GetValue(base + i);
     //set new-graph indices to Parent Graph nodes, other indices are set '-1'
     dest->SetValue(v, i);
     //this actually sets the 'data of Parent-Graph' in subgraph (undir)
     data->CopyData(this->Graph->GetVertexData(), v, i);
     }

  for (vtkIdType e = 0; e < this->Graph->GetNumberOfEdges(); e++)
     {
     u = this->Graph->GetSourceVertex(e);
     v = this->Graph->GetTargetVertex(e);
     u = dest->GetValue(u);
     v = dest->GetValue(v);

     if (u >= 0 & v >= 0 && undir->GetEdgeId(u, v) < 0)
       {
       builder->AddEdge(u, v);
       }
     }

  output->Initialize();
  output->ShallowCopy(builder->GetGraph());
  output->GetVertexData()->AddArray(nodeIds);//Adding JCN-ids in SubGraph
  output->Squeeze();
  output->Register(this);
  return output;	
}
//------------------------------------------------------------------------------------------
vtkSmartPointer<vtkIdTypeArray> vtkExtractReebGraph::GetCriticalNodes()
{
  if (!this->ReebGraph)
    {
    vtkErrorMacro("Input has not been processed yet, no ReebGraph available: First Extract the Reeb Graph..");
    return 0;
    }
  
  vtkDataSetAttributes *inPD   = this->ReebGraph->GetVertexData();
  vtkDataArray *inputjcnids = inPD->GetArray("jcnids");
  vtkDataArray *property = inPD->GetArray(this->FieldName);
  vtkDataArray *vertexids = inPD->GetArray("Vertex Ids");

  vtkIdType size = this->ReebGraph->GetNumberOfVertices();
  VTK_CREATE(vtkIdTypeArray, nodeIds);

  for (vtkIdType i = 0; i < size ; i++)
     {
     if (this->ReebGraph->GetDegree(i) == 2  )
       {
       //Find adjacent node values
       if ((this->ReebGraph->GetInDegree(i)==0 ||this->ReebGraph->GetOutDegree(i)==0))
         {
         nodeIds->InsertNextValue(inputjcnids->GetComponent(i,0));
         }
        }
      else //Not Regular Nodes 
        {
        nodeIds->InsertNextValue(inputjcnids->GetComponent(i,0));
        }
      }
   return nodeIds;
}

vtkReebGraph *vtkExtractReebGraph::GetReebGraph()
{
  if (!this->ReebGraph)
    {
    vtkErrorMacro("Input has not been processed, no graph available: First Extract the Reeb Graph..");
    return 0;
    }
  return ReebGraph;
}
//---------------ExtractReebGraph()-------------------------------------------------------
int vtkExtractReebGraph::ExtractReebGraph(vtkGraph *inputGraph, char *FieldName)
{
  //this->FieldName=FieldName;
  this->FieldName = new char[strlen(FieldName) + 1];
  strcpy(this->FieldName, FieldName);
  vtkDataSetAttributes *inPD   = inputGraph->GetVertexData();
  vtkDataArray *property = inPD->GetArray(FieldName); 

  if (!property)
    {
    vtkErrorMacro("Field not found in input.");
    return 0;
    }

  vtkIdType nrVertices = inputGraph->GetNumberOfVertices();
  //copying input graph JCN-ids
  vtkDataArray *inputjcnids;
  inputjcnids = inPD->GetArray("jcnids");

  //Setting JCN-ids
  VTK_CREATE(vtkIdTypeArray, nodeIds);
  nodeIds->SetName("jcnids");

  vtkSmartPointer<vtkDataArray> propCopy = property->NewInstance();
  propCopy->SetName(property->GetName());
  propCopy->SetNumberOfTuples(property->GetNumberOfTuples());

  // - union-find
  CreateUnionFindStructure(inputGraph, FieldName);

  // Generate output CT.   Map equivalence classes to nodes
  // in new graph, then translate edges in input graph into
  // edges between classes.
  VTK_CREATE(vtkMutableGraphHelper, builder); 
  //This is helpful in filters for (re)constructing graphs which may be either directed or undirected.
  VTK_CREATE(vtkMutableUndirectedGraph, undir);
  builder->SetGraph(undir);
  VTK_CREATE(vtkIdTypeArray, destination);
  destination->SetNumberOfValues(nrVertices);

  this->NrSupernodes = 0;	
  this->Size->SetNumberOfValues(nrVertices);  //Too large, but best we can do.

  //build tree nodes for each connected component
  for (vtkIdType i = 0; i < nrVertices; i++)
     {
     //This is equal to the number of connected components or the number of critical-nodes in CT
     if (this->UF->GetValue(i) < 0) 
       {
       destination->SetValue(i, this->NrSupernodes);
       propCopy->SetComponent(this->NrSupernodes, 0, property->GetComponent(i, 0));
       this->Size->SetValue(this->NrSupernodes, 0);
       builder->AddVertex();
       this->NrSupernodes++;
       //Copying Original JCN-id for ReebGraph-Nodes; 
       //Will be useful for the deepest level of MDRG;
        nodeIds->InsertNextValue(inputjcnids->GetComponent(i, 0));
       }
      }

  this->Start->SetNumberOfValues(this->NrSupernodes);
  this->Class->SetNumberOfValues(nrVertices);
  this->JCNIndex->SetNumberOfValues(nrVertices);

  vtkIdType root_i;
  for (vtkIdType i = 0; i < nrVertices; i++)
     {
     root_i = destination->GetValue(this->Root(i));
     destination->SetValue(i, root_i);		
     this->Size->SetValue(root_i, this->Size->GetValue(root_i) + 1);
     }
     
  vtkIdType u, v;
  // Copy across edges.
  for (vtkIdType e = 0; e < inputGraph->GetNumberOfEdges(); e++)
     {
     u = inputGraph->GetSourceVertex(e);
     v = inputGraph->GetTargetVertex(e);
     u = destination->GetValue(u);
     v = destination->GetValue(v);
     if (u != v && undir->GetEdgeId(u, v) < 0)
       {
       builder->AddEdge(u, v);
       }
     }
  VTK_CREATE(vtkIdTypeArray, base);
  base->SetNumberOfValues(this->NrSupernodes);

  this->Start->SetValue(0, 0);		
  base->SetValue(0, 0);
  // Compute absolute offset in list of field values for
  // each class	
  for (vtkIdType i = 1; i < this->NrSupernodes; i++)
     {
     this->Start->SetValue(i, this->Start->GetValue(i-1) + this->Size->GetValue(i-1));
     base->SetValue(i, this->Start->GetValue(i));
     }

  vtkIdType offset;
  // Store data for each node at offset.
  for (vtkIdType i = 0; i < nrVertices; i++)
     {
     root_i = destination->GetValue(i);
     offset = base->GetValue(root_i);
     this->Class->SetValue(offset, i);//Changing too next one?
     this->JCNIndex->SetValue(offset, inputjcnids->GetComponent(i, 0));
     base->SetValue(root_i, offset+1);
     }
  undir->ShallowCopy(builder->GetGraph()); 
  undir->GetVertexData()->AddArray(propCopy);

  //Creating ReebGraph Structure from this 'outputGraph'
  VTK_CREATE(vtkMutableDirectedGraph, dirGraph);
  vtkDataArray *outData = undir->GetVertexData()->GetArray(FieldName); 

  for (vtkIdType i = 0; i < undir->GetNumberOfVertices(); i++)
     {
     dirGraph->AddVertex();
     }

  for (vtkIdType i = 0; i < undir->GetNumberOfEdges(); i++)
     {
     vtkIdType u=undir->GetSourceVertex(i);
     vtkIdType v=undir->GetTargetVertex(i);
     double val_u=outData->GetComponent(u, 0);
     double val_v=outData->GetComponent(v, 0);
     if (val_u>val_v)
       {
       dirGraph->AddEdge(v, u);
       }
     else
       {
       dirGraph->AddEdge(u, v);
       }
     }
  dirGraph->GetVertexData()->AddArray(propCopy);
  dirGraph->GetVertexData()->AddArray(nodeIds);//Adding JCN-ids 
  VTK_CREATE(vtkReebGraph, outputReebGraph);
  bool success = outputReebGraph->CheckedShallowCopy(dirGraph);
  
  this->Graph = inputGraph;
  this->ReebGraph = outputReebGraph;
  return success;
}
//----------------------------------------------------------------------------

int vtkExtractReebGraph::RequestData(
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

  if (!this->FieldName)
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
     nodeIds->SetValue(i, i);//copy the origical JCN-ids!	
     }
  inputGraph->GetVertexData()->AddArray(nodeIds);//Adding JCN-ids in Input Graph
  ExtractReebGraph(inputGraph, this->FieldName);

  //Designing OutputGraph
  vtkDataSetAttributes *inPD   = inputGraph->GetVertexData(), 
                       *outPD  = outputGraph->GetVertexData();
  outputGraph->ShallowCopy(inputGraph);
  outPD->CopyAllocate(inPD, nrVertices, 10);	
  for (vtkIdType i = 0; i < nrVertices; i++)
     {
     outPD->CopyData(inputGraph->GetVertexData(), i, i);
     }

  outputGraph->Squeeze();	
  return 1;
}
//-------------------------------------------------------------------------------------
 
