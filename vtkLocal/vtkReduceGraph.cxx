/*=========================================================================

 *	File: vtkReduceGraph.cxx
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
#include "vtkReduceGraph.h"
#include "vtkTreeBFSIterator.h"
#include "vtkVoidArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkInformationVector.h"
#include "vtkExecutive.h"

vtkStandardNewMacro(vtkReduceGraph);

#define VTK_CREATE(type,name) \
   vtkSmartPointer<type> name = vtkSmartPointer<type>::New()


//-----------------------------------------------------------------------
vtkReduceGraph::vtkReduceGraph()
{
  this->SetNumberOfOutputPorts(2);
  this->SetNumberOfInputPorts(2);
  
  vtkUndirectedGraph *output1 = vtkUndirectedGraph::New();
  this->GetExecutive()->SetOutputData(0, output1);
  output1->Delete();

  vtkUnstructuredGrid *output2 = vtkUnstructuredGrid::New();
  this->GetExecutive()->SetOutputData(1, output2);
  output2->Delete();
}

vtkReduceGraph::~vtkReduceGraph()
{
}



void vtkReduceGraph::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

int vtkReduceGraph::FillOutputPortInformation(int port, vtkInformation* info)
{
  // Port 1 = optional unstructured grid of slabs.
  if (port == 0)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUndirectedGraph");
    return 1;
    }
  if (port == 1)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
    }
  return 0;
}

int vtkReduceGraph::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0)
    {
     info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkGraph");
     return 1;
    }
  
  if (port == 1)
    {
     info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
     return 1;
    }
  return 0;
}

vtkDataObject* vtkReduceGraph::GetOutput(int index)
{
  vtkDataObject *obj;
  if (index == 0)
    {
    obj = vtkUndirectedGraph::SafeDownCast(this->GetOutputDataObject(0));
    return obj;
    }
  else if (index == 1)
    {
    obj = vtkUnstructuredGrid::SafeDownCast(this->GetOutputDataObject(1));
    return obj;
    }
  return NULL;
}

int vtkReduceGraph::RequestDataObject(
  vtkInformation* request,
  vtkInformationVector** inputVector ,
  vtkInformationVector* outputVector)
{
  vtkInformation* info0 = outputVector->GetInformationObject(0);
  vtkInformation* info1 = outputVector->GetInformationObject(1);

  vtkUndirectedGraph *output = vtkUndirectedGraph::SafeDownCast(
        info0->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *outputGrid = vtkUnstructuredGrid::SafeDownCast(
        info1->Get(vtkDataObject::DATA_OBJECT()));

  if (!output)
    {
    output = vtkUndirectedGraph::New();
    info0->Set(vtkDataObject::DATA_OBJECT(), output);
    output->Delete();
    }

  if (!outputGrid)
    {
    outputGrid = vtkUnstructuredGrid::New();
    info1->Set(vtkDataObject::DATA_OBJECT(), outputGrid);
    outputGrid->Delete();
    }
    
  return 1;
}

//----------------------------------------------------------------------------
int vtkReduceGraph::RequestData(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector, 
  vtkInformationVector* outputVector)
{
  vtkInformation *inInfo1 = inputVector[0]->GetInformationObject(0);
  vtkInformation *inInfo2 = inputVector[1]->GetInformationObject(0);
  
  vtkInformation *outInfo1 = outputVector->GetInformationObject(0);
  vtkInformation *outInfo2 = outputVector->GetInformationObject(1);
  
    // Output graph data
  vtkGraph *inputGraph = vtkGraph::SafeDownCast(
         inInfo1->Get(vtkDataObject::DATA_OBJECT()));
         
  vtkUnstructuredGrid *inputGrid = vtkUnstructuredGrid::SafeDownCast(
         inInfo2->Get(vtkDataObject::DATA_OBJECT()));
  
    // Output graph data
  vtkGraph *outputGraph = vtkGraph::SafeDownCast(
         outInfo1->Get(vtkDataObject::DATA_OBJECT()));
         
  vtkUnstructuredGrid *outputGrid = vtkUnstructuredGrid::SafeDownCast(
         outInfo2->Get(vtkDataObject::DATA_OBJECT()));
         
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
 
   if (!inputGrid)
     {
     vtkErrorMacro("A vtkUnstructuredGrid object is required.");
     return 0;
     }

  vtkIdType nrCells = inputGrid->GetNumberOfCells();
  if (!nrCells)
    {
    vtkErrorMacro("Empty unstructured grid, no cells.");
    }

  vtkSmartPointer<vtkIdTypeArray> idArray = vtkSmartPointer<vtkIdTypeArray>::New();
  idArray->SetName("JcnId");
  for (int i = 0; i < inputGraph->GetNumberOfVertices(); i++)
    {
    idArray->InsertNextValue(i);
    }
    
  inputGraph->GetVertexData()->AddArray(idArray);
  
  VTK_CREATE(vtkBoostConnectedComponents, comp);
  comp->SetInputData(inputGraph);
  comp->Update();
  
  vtkSmartPointer<vtkGraph> compGraph = comp->GetOutput();
  int compNr = compGraph->GetVertexData()->GetArray("component")->GetRange()[1]+1;

  int GetComponentNr[compNr];
  for (int i = 0; i < compNr; i++) GetComponentNr[i] = 0;
  
  for (int i = 0; i < compGraph->GetVertexData()->GetArray("component")->GetNumberOfTuples();
       i++)
     {
     int nr = compGraph->GetVertexData()->GetArray("component")->GetComponent(i,0);
     GetComponentNr[nr]++;
     }
 
  int maxCompNr = 0;
  for (int i = 1; i < compNr; i++)
     {
     if (GetComponentNr[i] > GetComponentNr[maxCompNr]) maxCompNr = i;
     }

  vtkSmartPointer<vtkIdTypeArray> idList = vtkSmartPointer<vtkIdTypeArray>::New();
  for (int i = 0; i < compGraph->GetNumberOfVertices(); i++)
     {
     int nr = compGraph->GetVertexData()->GetArray("component")->GetComponent(i,0);
     if (nr == maxCompNr) idList->InsertNextValue(i);
     }
  
  vtkSmartPointer<vtkSelectionNode> selectionNode = vtkSmartPointer<vtkSelectionNode>::New();	 
  selectionNode->SetFieldType(vtkSelectionNode::VERTEX);
  selectionNode->SetContentType(vtkSelectionNode::INDICES);
  vtkSmartPointer<vtkSelection> selection = vtkSmartPointer<vtkSelection>::New();
  vtkSmartPointer<vtkExtractSelectedGraph> extractSelection = vtkSmartPointer<vtkExtractSelectedGraph>::New();
  
  selectionNode->SetSelectionList(idList);  
  selection->AddNode(selectionNode);
  extractSelection->SetInputData(0, inputGraph);
  extractSelection->SetInputData(1, selection);
  extractSelection->Update();
  outputGraph->ShallowCopy(extractSelection->GetOutput());
   
  vtkSmartPointer<vtkIdTypeArray> jcnIdList = vtkSmartPointer<vtkIdTypeArray>::New();
  for (int i = 0; i < outputGraph->GetNumberOfVertices(); i++)
     {
     vtkIdType id = outputGraph->GetVertexData()->GetArray("JcnId")->GetComponent(i,0);
     jcnIdList->InsertNextValue(id);
     }
  
  std::cout << "OutGraph size: " << outputGraph->GetNumberOfVertices() << std::endl;
  std::cout << "jcnIdList size: "<< jcnIdList->GetNumberOfTuples() << std::endl;
  
  
  selectionNode->SetFieldType(vtkSelectionNode::CELL);
  selectionNode->SetContentType(vtkSelectionNode::INDICES);
  selectionNode->SetSelectionList(jcnIdList);
  selection->AddNode(selectionNode);
  vtkSmartPointer<vtkExtractSelection> extractSelectionGrid = vtkSmartPointer<vtkExtractSelection>::New();
  extractSelectionGrid->SetInputData(0, inputGrid);
  extractSelectionGrid->SetInputData(1, selection);
  extractSelectionGrid->Update();

  outputGrid->ShallowCopy(extractSelectionGrid->GetOutput());
  outputGrid->Squeeze();
  
  return 1;
}
//--------------------------------------------------------------------------------------
