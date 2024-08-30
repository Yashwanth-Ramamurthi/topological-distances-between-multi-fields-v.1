/*=========================================================================

 *	File: vtkSimilarityBetweenMultiResolutionReebSpaces.h
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/

// .NAME vtkSimilarityBetweenMultiResolutionReebSpaces - Computes the similarity between multi-resolution Reeb spaces of two multi-fields
//Input: Joint Countour Nets (finest resolution Reeb spaces) of two multi-fields
//Output: Similarity between the multi-resolution Reeb spaces

//Reference:
//[1] Yashwanth Ramamurthi and Amit Chattopadhyay. "A Topological Similarity Measure Between Multi-Resolution Reeb Spaces," in IEEE Transactions on Visualization and Computer Graphics, vol. 28, no. 12, pp. 4360-4374, 1 Dec. 2022

#ifndef __vtkSimilarityBetweenMultiResolutionReebSpaces_h
#define __vtkSimilarityBetweenMultiResolutionReebSpaces_h
#include "vtkMultiResolutionReebSpace.h"
#include "vtkDataSetAttributes.h"
#include "vtkIdTypeArray.h"
#include "vtkBitArray.h"
#include "vtkFloatArray.h"
#include "vtkIntArray.h"
#include "vtkInformation.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkMutableUndirectedGraph.h"
#include "vtkMutableGraphHelper.h"
#include "vtkAdjacentVertexIterator.h"
#include "vtkGraphEdge.h"
#include "vtkObjectFactory.h"
#include "vtkBoostConnectedComponents.h"
#include "vtkNew.h"
#include <queue>
#include <stack>
#include <algorithm>
#include <set>
#include <vector>
#include "vtkVoidArray.h"
#include "vtkTreeBFSIterator.h"
#include "vtkDoubleArray.h"
#include <vtkExtractReebGraph.h>
#include "vtkMultiResolutionReebSpace.h"
#include <vtkVariant.h>
#include <vtkVariantArray.h>
using namespace std;

class VTK_META_EXPORT vtkSimilarityBetweenMultiResolutionReebSpaces : public vtkGraphAlgorithm
{
public:
clock_t start, mid, end;
int fileName1, fileName2;
	//Priority queue consisting the nodes in MRS1 where the priority of a node is its volume (size)
	priority_queue<pair<double, int> > queueMrs1Nodes;
	
	//outer vector for resolution, inner vector for matched node pairs between JCNs of a particular resolution
	vector < vector<pair <int,int> > > MPAIRS;

	double totalSizeMrs1, totalSizeMrs2;
	double weightVolume;
	double weightRangeMeasure;
	double weightBettiZero;
	double weightDegree;
	vtkMultiResolutionReebSpace mrs1, mrs2;
	vtkMultiDimensionalReebGraph *mdrg1, *mdrg2;
	map<int, vector<pair<vtkGraph*, int> > > mdrgNodeMap1, mdrgNodeMap2;
	map<pair<vtkGraph*, int>, vector<int> > mlistMap1, mlistMap2;
	vector<bool> IsMatchedMrs1;
	vector<bool> IsMatchedMrs2;
	int numberOfResolutions;
	vector<char*> Fields;
	std::vector<double> minimumValues;
	std::vector<double> maximumValues;
	std::vector<double> slabwidths;

	static vtkSimilarityBetweenMultiResolutionReebSpaces* New();

	vtkTypeMacro(vtkSimilarityBetweenMultiResolutionReebSpaces, vtkGraphAlgorithm);

	vtkSimilarityBetweenMultiResolutionReebSpaces();

	~vtkSimilarityBetweenMultiResolutionReebSpaces();
	
	//Set weights for the four attributes required for computing the similarity measure (See Equation (3) in [1])
	void SetAttributes(double weightVolume, double weightRangeMeasure, double weightBettiZero, double weightDegree);

	//Create maps between the nodes of the JCN and the MDRG, and create an empty mlist for each of the MDRG nodes
	void createMdrgNodeMap(int graphId, vector<pair<vtkGraph*, int> > nodes, int mrsId);

	//Initialize the data-structures for computing the similarity between two multi-resolution Reeb spaces
	void Initialize(int resolution);

	//Compute the similarity between two nodes based on the the volumes (sizes) of the corresponding quantized fiber-components (See Sections 4.3 and 5.4 in [1])
	double SimilarityBetweenVolumes(int nodeId1, int nodeId2, int resolution);

	//Compute the similarity between two nodes based on the their values in the range (See Sections 4.3 and 5.4 in [1])
	double SimilarityBetweenRanges(int nodeId1, int nodeId2, int resolution);

	//Compute the similarity between two nodes based on the number of components of the corresponding joint level sets (See Sections 4.3 and 5.4 in [1])
	double SimilarityBetweenNumberOfComponents(int nodeId1, int nodeId2, int resolution);

	//Compute the similarity between two nodes based on the their degrees (See Sections 4.3 and 5.4 in [1])
	double SimilarityBetweenDegrees(int nodeId1, int nodeId2, int resolution);

	//Compute the similarity between two nodes (See Equation (3) in [1])
	double SimilarityBetweenNodes(int nodeId1, int nodeId2, int resolution);
	
	//Compute the similarity between the MPAIRS of two multi-resolution Reeb Spaces
	double SimilarityBetweenMPAIRS();

	//Check if the parents of the nodes to be matched form an MPAIR (See Section 5.1 in [1])
	bool CheckParentMatching(int parentID1, int parentID2);

	//Check whether the mlists of the nodes to be matched are the same (topological consistency) - See Section 5.1 in [1]
	bool CheckMLIST(vector<int> mlist1, vector<int> mlist2);

	//Check the matching criteria (See Section 5.1 in [1])
	bool CheckMatchingCriteria(int nodeId1, int nodeId2, int resolution);

	//Find candidate nodes satisfying the matching criteria (See Section 5.1 in [1])
	vector<int> FindCandidateMatchingNodes(int nodeId, int resolution);

	//Propogate the label along a branch in the monotonically increasing or decreasing direction (See Section 5.2 in [1])
	void PropagateLabel(pair<vtkGraph*, int> node, int matchingLabel, bool isDecreasing, int reebSpace);

	//Finding the MPAIRs (matching pairs of nodes) between JCNs of a particular resolution (see Algorithm 2 in [1])
	void FindMatchingPairs(int resolution);

	//Add the parameters for computing the two multi-resolution Reeb spaces and the similarity between them
	void AddParameters(char *fieldName, double minimumValue,double maximumValue,double slabwidth);

	//Compute the total size of the fiber-components corresponding to all the nodes in each of the multi-resolution Reeb spaces
	void ComputeTotalSize();

	//Compute the similarity between two multi-resolution Reeb spaces with finest resolution JCNs "JCN1" and "JCN2", respectively.
	double ComputeSimilarity(vtkGraph * JCN1, vtkGraph * JCN2, int numberOfResolutions);

};
#endif
