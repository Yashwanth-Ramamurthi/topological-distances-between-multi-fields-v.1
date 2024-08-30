/*=========================================================================

 *	File: vtkMultiResolutionReebSpace.h
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/

// .NAME vtkMultiResolutionReebSpace - Computes the multi-resolution Reeb space of a multi-field f = (f1,f2,...,fn)
//Input: Joint Countour Net (finest resolution Reeb space) corresponding to a multi-field
//Output: Multi-resolution Reeb space

//Reference:
//[1] Yashwanth Ramamurthi and Amit Chattopadhyay. "A Topological Similarity Measure Between Multi-Resolution Reeb Spaces," in IEEE Transactions on Visualization and Computer Graphics, vol. 28, no. 12, pp. 4360-4374, 1 Dec. 2022

#ifndef __vtkMultiResolutionReebSpace_h
#define __vtkMultiResolutionReebSpace_h
#include "vtkMETAModule.h"
#include "vtkBitArray.h"
#include "vtkSmartPointer.h"
#include "vtkGraphAlgorithm.h"
#include <vtkMutableGraphHelper.h>
#include <vtkExtractReebGraph.h>
#include <vector>
#include "vtkMultiDimensionalReebGraph.h"
#include "vtkSelectionNode.h"
#include "vtkSelection.h"
#include "vtkExtractSelection.h"
#include "vtkExtractSelectedGraph.h"
#include "vtkAdjacentVertexIterator.h"
#include "vtkTimerLog.h"
#include<map>
using namespace std;
class vtkIdTypeArray;
class VTK_META_EXPORT vtkMultiResolutionReebSpace : public vtkGraphAlgorithm
{
public:

static vtkMultiResolutionReebSpace* New();
void PrintSelf(ostream& os, vtkIndent indent);
vtkMultiResolutionReebSpace();
~vtkMultiResolutionReebSpace();

//Find the nodes in the finer resolution JCN which will be merged to a single node in the coarser resolution JCN
//Nodes in "component" (returned by this procedure) will be merged to a single node in the coarser resolution JCN
vector<int> UnionFind(vtkGraph *jcn, int nodeId, vector<bool> &visited,vector<int> component, vector<vector<pair<double,double> > > rangePairs);

//Obtain the id of the component containing the node with id "nodeId"
int GetComponentNumber(int nodeId, vector<vector<int> > components);

//Compute a coarser JCN from a fine JCN. Nodes in the finer JCN will be merged to obtain the coarser JCN
vtkSmartPointer<vtkUndirectedGraph> ComputeCoarseJCN(vtkGraph *fineJCN,vector<double> minimumValue, vector<double> maximumValue, vector<double> slabwidth);

//Construct the multi-resolution Reeb Space by taking the JCN in the finest resolution as input
void ConstructMrs(vtkGraph* finestResolutionJCNGraph, vector <char*> fields,vector<double> minimumValue, vector<double> maximumValue, vector<double> slabwidth, int numberOfResolutions);

//Stores the names of the component scalar fields of the multi-fields (eg., f1,f2,..,fn).
std::vector<char*> Fields;

//Stores the quantized Reeb spaces (JCNs) in multiple resolutions
std::vector< vtkSmartPointer < vtkUndirectedGraph> > jcnGraphs;

//Stores the list of nodes corresponding every combination of range values of the fields for JCNs of multiple resolutions
vector< map<vector<double>, vector<int> > > nodeRangeMaps;

//Stores the minimum ranges for each field of each JCN
vector<vector<double> > minimumRanges;

//Stores the maximum ranges for each field of each JCN 
vector<vector<double> > maximumRanges;

//Number of resolutions of the multi-resolution Reeb space
int numberOfResolutions;

private:
vtkMultiResolutionReebSpace(const vtkMultiResolutionReebSpace&);   // Not implemented
void operator=(const vtkMultiResolutionReebSpace&);   // Not implemented
};

#endif
