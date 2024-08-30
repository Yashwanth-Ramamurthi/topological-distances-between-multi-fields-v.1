/*=========================================================================

 *	File: vtkExtractReebGraph.h
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
// .NAME vtkExtractReebGraph - "Extracts" the Reeb Graph from a subgraph of the JCN.
// .SECTION Description
// vtkExtractReebGraph extracts a Reeb graph from a subgraph of the input JCN. 
// This filter takes the computed JCN graph as input. To initialise the structure, 
// names of each scalar field has to be given. 
// Conceptually, we first create a Union-Find structure and compute the connected 
// components of field f. Then we create Reeb graph nodes corresponding to the connected 
// components and order them according to the value of f. Finally, create a "directed" (from
// lower f-value to higher f-value) edge between two nodes if they are adjacent components 
// in the JCN.


#ifndef __vtkExtractReebGraph_h
#define __vtkExtractReebGraph_h

#include "vtkMETAModule.h"
#include "vtkSmartPointer.h"
#include "vtkGraphAlgorithm.h"
#include <vtkMutableGraphHelper.h>
#include <vtkReebGraph.h>

class vtkIdTypeArray;

class VTK_META_EXPORT vtkExtractReebGraph : public vtkGraphAlgorithm
{
public:
  static vtkExtractReebGraph* New();
  vtkTypeMacro(vtkExtractReebGraph,vtkGraphAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
 
  // Description: 
  // Get/Set the data field name
  vtkSetStringMacro(FieldName);
  vtkGetStringMacro(FieldName);
  
  // Description:
  // Return the subgraph of a group of selected graph Ids. 
  vtkGraph *GetSubGraph(vtkIdType groupId);
  
  // Description: 
  // Return the size of the subgraph of selected graph ids. 
  vtkIdType GetSizeOfSubGraph(vtkIdType groupId);
  
  // Description: 
  // Return the JCN ids of the selected subgraph. 
  vtkSmartPointer<vtkIdTypeArray> GetSubGraphJCNIds(vtkIdType groupId);
  
  // Description:
  // ReebGraph-functions
  // Creates a Reeb graph from an input subgraph of the input JCN and a field-name
  int ExtractReebGraph(vtkGraph *inputGraph, char *FieldName);
  
  // Description:
  // Returns the created Reeb graph
  vtkReebGraph *GetReebGraph();
  
  // Description:
  // Returns critical nodes of Reeb Graph
  vtkSmartPointer<vtkIdTypeArray> GetCriticalNodes();

protected:
  vtkExtractReebGraph();
  ~vtkExtractReebGraph();

  int FillInputPortInformation(int port, vtkInformation* info);
  
  int RequestData(
    vtkInformation*, 
    vtkInformationVector**, 
    vtkInformationVector*);
  
  // Find the root of a value stored in the union-find
  // structure, performing path compression as we go.
  vtkIdType Root(vtkIdType x);
  
  // Find if x and y are in the same component
  bool Find(vtkIdType x, vtkIdType y);
  
  //Merge components containing i and nb, if they are not in a same component
  void Union(vtkIdType p, vtkIdType q);
  void CreateUnionFindStructure(vtkGraph *inputGraph, char *FieldName);
  
  // List of names for the scalar fields
  char *FieldName;

  // Input Graph
  vtkSmartPointer<vtkGraph> Graph;

  // Output Graph
  vtkSmartPointer<vtkReebGraph> ReebGraph;

  // union-find structure
  vtkSmartPointer<vtkIdTypeArray> UF;	

  // data for Extracted Sub-Graphs
  vtkSmartPointer<vtkIdTypeArray> Class;

  //Parent graph node ids; Usage: pid(i) = this->Class->GetValue(base + i)
  vtkSmartPointer<vtkIdTypeArray> JCNIndex;

  //Parent graph node ids; Usage: pid(i) = this->Class->GetValue(base + i)
  //base ids of Sub-Graphs; Usage: base = this->Start->GetValue(groupId) 
  vtkSmartPointer<vtkIdTypeArray> Start;
  
  //Sizes of Extracted Sub-Graphs; Usage: size = this->Size->GetValue(groupId)
  vtkSmartPointer<vtkIdTypeArray> Size;
  vtkIdType NrSupernodes;

private:
  vtkExtractReebGraph(const vtkExtractReebGraph&); // Not implemented
  void operator=(const vtkExtractReebGraph&);   // Not implemented
};

#endif

