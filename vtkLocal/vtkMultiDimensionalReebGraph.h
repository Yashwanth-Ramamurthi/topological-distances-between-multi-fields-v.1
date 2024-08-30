/*=========================================================================

 *	File: vtkMultiDimensinoalReebGraph.h
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/

// .NAME vtkMultiDimensionalReebGraph - Creates a Multi-Dimensional Reeb Graph from the JCN.
// .SECTION Description
// Conceptually, it calls the vtkExtractReebGraph, recursively, to create the 
// Multi-Dimensional Reeb graph structure. It takes one input - the computed JCN graph
// structure. The output is the same JCN graph, but with an additional array indicating 
// if each node is part of Jacobi Sets or not. 
// The main algorithm is based on the following publication : 
//
// Amit Chattopadhyay, Hamish Carr, David Duke and Zhao Geng, Extracting Jacobi Structures 
// in Reeb Spaces. The Eurographics Conference on Visualization, Eurovis ShortPaper 2014, 
// 9-13 June, Swansea, UK


#ifndef __vtkMultiDimensionalReebGraph_h
#define __vtkMultiDimensionalReebGraph_h

#include "vtkExtractReebGraph.h"
#include "vtkGraphAlgorithm.h"
#include "vtkMETAModule.h"
#include "vtkMutableGraphHelper.h"
#include "vtkSmartPointer.h"
#include "vtkTree.h"
#include <vector>

class vtkIdTypeArray;

class VTK_META_EXPORT vtkMultiDimensionalReebGraph : public vtkGraphAlgorithm
{
public:
  static vtkMultiDimensionalReebGraph* New();
  vtkTypeMacro(vtkMultiDimensionalReebGraph, vtkGraphAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  // Description:
  // Initialise the names for each field.  
  void AddField(char *nm);
  // Description:  
  // MDRG-functions
  // Creates a MultiDimensionalReebGraph structure and stores in a vtkTree structure
  void CreateMultiDimensionalReebGraph(vtkGraph *inputGraph);
  // Description:
  // Returns the MultiDimensionalReebGraph tree-structure
  vtkTree *GetMultiDimensionalReebGraph();
  // Description:
  // Get a tree node (a Reeb Graph) from the MultiDimensionalReebGraph tree-structure
  vtkSmartPointer<vtkExtractReebGraph> GetTreeNode(vtkIdType i);
  // Description:
  // This function extracts the Jacobi Nodes from the MultiDimensionalReebGraph tree-structure.
  // These are the critical nodes of the Reeb Graphs at the deepest level of the MultiDimensionalReebGraph
  vtkSmartPointer<vtkIdTypeArray> ExtractJacobiNodes();

protected:
  vtkMultiDimensionalReebGraph();
  ~vtkMultiDimensionalReebGraph();

  int FillInputPortInformation(int port, vtkInformation* info);
  
  int RequestData(
    vtkInformation*, 
    vtkInformationVector**, 
    vtkInformationVector*);
  // Field names
  vtkstd::vector<char*> Fields;
  // Jacobi Structure
  vtkSmartPointer<vtkTree> MDRG;
  // Computed nested Reeb Graph structures. 
  std::vector<vtkSmartPointer<vtkExtractReebGraph> > ReebGraphList;


private:
  vtkMultiDimensionalReebGraph(const vtkMultiDimensionalReebGraph&); // Not implemented
  void operator=(const vtkMultiDimensionalReebGraph&);   // Not implemented
};

#endif

