/*=========================================================================

 *	File: vtkReduceGraph.h
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
// .NAME vtkReduceGraph - Given a disconnected VTK graph, this filter extracts the 
// largest connected subgraph . 
// .SECTION Description
//  The input is a vtkGraph object and the output is the largest connected subgraph of 
//  the input graph. 
#ifndef __vtkReduceGraph_h
#define __vtkReduceGraph_h
#include "vtkAdjacentVertexIterator.h"
#include "vtkBitArray.h"
#include "vtkDijkstraGraphGeodesicPath.h"
#include "vtkExtractSelection.h"
#include "vtkExtractSelectedGraph.h"
#include "vtkGraphAlgorithm.h"
#include "vtkGraphToPolyData.h"
#include "vtkIdList.h"
#include "vtkMETAModule.h"
#include "vtkMutableGraphHelper.h"
#include "vtkMultiDimensionalReebGraph.h"
#include "vtkSmartPointer.h"
#include "vtkSelectionNode.h"
#include "vtkSelection.h"
#include <ctime>
#include <vector>
class vtkIdTypeArray;

class VTK_META_EXPORT vtkReduceGraph : public vtkGraphAlgorithm
{
public:
  static vtkReduceGraph* New();
  vtkTypeMacro(vtkReduceGraph, vtkGraphAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  vtkDataObject* GetOutput(int index);
protected:
  vtkReduceGraph();
  ~vtkReduceGraph();
  
  int FillInputPortInformation(int port, vtkInformation* info);
  int FillOutputPortInformation(int port, vtkInformation* info);
  int RequestDataObject(vtkInformation* request,vtkInformationVector** inputVector,
                        vtkInformationVector* outputVector);
  
  int RequestData(
    vtkInformation*, 
    vtkInformationVector**, 
    vtkInformationVector*);
  
  private:
  vtkReduceGraph(const vtkReduceGraph&);   // Not implemented
  void operator=(const vtkReduceGraph&);   // Not implemented
};

#endif

