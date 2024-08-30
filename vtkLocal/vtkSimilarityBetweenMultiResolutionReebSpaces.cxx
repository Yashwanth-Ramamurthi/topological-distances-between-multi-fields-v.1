/*=========================================================================

 *	File: vtkSimilarityBetweenMultiResolutionReebSpaces.cxx
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
//Reference:
//[1] Yashwanth Ramamurthi and Amit Chattopadhyay. "A Topological Similarity Measure Between Multi-Resolution Reeb Spaces," in IEEE Transactions on Visualization and Computer Graphics, vol. 28, no. 12, pp. 4360-4374, 1 Dec. 2022

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
#include <cmath>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include "vtkSimilarityBetweenMultiResolutionReebSpaces.h"
#include <vtkInEdgeIterator.h>
#include <vtkOutEdgeIterator.h>
#include <math.h>
#include<bits/stdc++.h>
#include<stdlib.h>
#include<stdio.h>
#include<fstream>
#include<sstream>
using namespace std;

vtkStandardNewMacro(vtkSimilarityBetweenMultiResolutionReebSpaces);

#define VTK_CREATE(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

vtkSimilarityBetweenMultiResolutionReebSpaces::vtkSimilarityBetweenMultiResolutionReebSpaces(){
    this->Fields = std::vector<char*>();
}

vtkSimilarityBetweenMultiResolutionReebSpaces::~vtkSimilarityBetweenMultiResolutionReebSpaces(){
}

//Set weights for the four attributes required for computing the similarity measure (See Equation (3) in [1])
void vtkSimilarityBetweenMultiResolutionReebSpaces::SetAttributes(double weightVolume, double weightRangeMeasure, double weightBettiZero, double weightDegree){
	this->weightVolume = weightVolume;
	this->weightRangeMeasure = weightRangeMeasure;
	this->weightBettiZero = weightBettiZero;
	this->weightDegree = weightDegree;
}

//Create maps between the nodes of the JCN and the MDRG, and create an empty mlist for each of the MDRG nodes
void vtkSimilarityBetweenMultiResolutionReebSpaces::createMdrgNodeMap(int graphId, vector<pair<vtkGraph*, int> > nodes, int mrsId){
	vtkGraph *reebGraph;
	vtkMultiDimensionalReebGraph *mdrg;
	if(mrsId == 1){
		mdrg = mdrg1;
	}
	else{
		mdrg = mdrg2;
	}
mdrg->GetTreeNode(graphId);
	reebGraph = mdrg->GetTreeNode(graphId)->GetReebGraph();
	if(mdrg->GetMultiDimensionalReebGraph()->IsLeaf(graphId)){
		int nrVertices = reebGraph->GetNumberOfVertices();
		vtkIdTypeArray *jcnIdArray = vtkIdTypeArray::SafeDownCast(reebGraph->GetVertexData()->GetArray("jcnids"));
		for(int i = 0; i < nrVertices; i++){
			int jcnId = jcnIdArray->GetComponent(i,0);
			vector<pair<vtkGraph*, int> > nodes1 = nodes;
			nodes1.push_back(make_pair(reebGraph, i));
			if(mrsId == 1){
				mdrgNodeMap1[jcnId] = nodes1;
				vector<int> mlist;
				mlistMap1[make_pair(reebGraph, i)] = mlist;
			}
			else{
				mdrgNodeMap2[jcnId] = nodes1;	
				vector<int> mlist;
				mlistMap2[make_pair(reebGraph, i)] = mlist;
			}
		}
	}
	else{
		vtkSmartPointer<vtkAdjacentVertexIterator> adj = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
		mdrg->GetMultiDimensionalReebGraph()->GetAdjacentVertices(graphId, adj);
		int vertexId = 0;
		while (adj->HasNext()){
			vector<pair<vtkGraph*, int> > nodes1 = nodes;
			nodes1.push_back(make_pair(reebGraph, vertexId));
			int adjacentVertex = adj->Next();
			createMdrgNodeMap(adjacentVertex, nodes1, mrsId);
			vertexId++;
		}
	}
}

//Initialize the data-structures for computing the similarity between two multi-resolution Reeb spaces
void vtkSimilarityBetweenMultiResolutionReebSpaces::Initialize(int resolution){
	mdrgNodeMap1.clear();
	mdrgNodeMap2.clear();
	IsMatchedMrs1.clear();
	IsMatchedMrs2.clear();
	mdrg1 = vtkMultiDimensionalReebGraph::New();
	for(int f = 0; f < this->Fields.size(); f++){
		mdrg1->AddField(this->Fields[f]);
	}
	mdrg1->CreateMultiDimensionalReebGraph(mrs1.jcnGraphs[resolution]);
	mdrg2 = vtkMultiDimensionalReebGraph::New();
	for(int f = 0; f < this->Fields.size(); f++){
		mdrg2->AddField(this->Fields[f]);
	}
	mdrg2->CreateMultiDimensionalReebGraph(mrs2.jcnGraphs[resolution]);
	vector<pair<vtkGraph*, int> > nodes1, nodes2;
	createMdrgNodeMap(0, nodes1, 1);
	createMdrgNodeMap(0, nodes2, 2);
	int nrVertices = mrs1.jcnGraphs[resolution]->GetNumberOfVertices();
	vtkIntArray *sizeArray = vtkIntArray::SafeDownCast(mrs1.jcnGraphs[resolution]->GetVertexData()->GetArray("size"));
	for(int i = 0; i < nrVertices; i++){
			IsMatchedMrs1.push_back(false);
	}
	if(resolution == 0){
		for(int i = 0; i < nrVertices; i++){
			queueMrs1Nodes.push(make_pair(sizeArray->GetComponent(i,0), -i));
		}
	}
	else{
		vtkSmartPointer<vtkVariantArray> children = vtkVariantArray::SafeDownCast(mrs1.jcnGraphs[resolution-1]->GetVertexData()->GetAbstractArray("children"));
		for(int i = 0; i < this->MPAIRS[resolution - 1].size(); i++){
			int parentNodeId = this->MPAIRS[resolution - 1][i].first;
			vtkIdTypeArray* childArray = vtkIdTypeArray::SafeDownCast(children->GetValue(parentNodeId).ToArray());
			int numberOfVertices = childArray->GetNumberOfTuples();
			for(int j = 0; j < numberOfVertices; j++){
				int nodeId = childArray->GetComponent(j,0);
				queueMrs1Nodes.push(make_pair(sizeArray->GetComponent(nodeId,0), -nodeId));
			}
		}
	}
	nrVertices = mrs2.jcnGraphs[resolution]->GetNumberOfVertices();
	for(int i = 0; i < nrVertices; i++){
		IsMatchedMrs2.push_back(false);
	}
}

//Compute the similarity between two nodes based on the the volumes (sizes) of the corresponding quantized fiber-components (See Sections 4.3 and 5.4 in [1])
double vtkSimilarityBetweenMultiResolutionReebSpaces::SimilarityBetweenVolumes(int nodeId1, int nodeId2, int resolution){
	double size1 =  (double)mrs1.jcnGraphs[resolution]->GetVertexData()->GetArray("size")->GetComponent(nodeId1,0)/(double)totalSizeMrs1;
	double size2 =  (double)mrs2.jcnGraphs[resolution]->GetVertexData()->GetArray("size")->GetComponent(nodeId2,0)/(double)totalSizeMrs2;
	return min(size1,size2)/max(size1,size2);
}

//Compute the similarity between two nodes based on the their values in the range (See Sections 4.3 and 5.4 in [1])
double vtkSimilarityBetweenMultiResolutionReebSpaces::SimilarityBetweenRanges(int nodeId1, int nodeId2, int resolution){
	double valueMrs1 = 1.0, valueMrs2 = 1.0;
	for(int f = 0; f < this->Fields.size(); f++){
		double slabwidth = this->slabwidths[f]*pow(2, numberOfResolutions - 1 - resolution);
		valueMrs1 = valueMrs1 * ( slabwidth/(mrs1.maximumRanges[resolution][f] + slabwidth - mrs1.minimumRanges[resolution][f]));
		valueMrs2 = valueMrs2 * ( slabwidth/(mrs2.maximumRanges[resolution][f] + slabwidth - mrs2.minimumRanges[resolution][f]));
	}
	return min(valueMrs1,valueMrs2)/max(valueMrs1,valueMrs2);
}

//Compute the similarity between two nodes based on the number of components of the corresponding joint level sets (See Sections 4.3 and 5.4 in [1])
double vtkSimilarityBetweenMultiResolutionReebSpaces::SimilarityBetweenNumberOfComponents(int nodeId1, int nodeId2, int resolution){
	vector<double> rangeVector;
	vtkDataSetAttributes *graphProperties = mrs1.jcnGraphs[resolution]->GetVertexData();
	for(int f = 0; f < this->Fields.size(); f++){
		rangeVector.push_back(graphProperties->GetArray(this->Fields[f])->GetComponent(nodeId1,0));
	}
	double numberComponents_mrs1 = mrs1.nodeRangeMaps[resolution][rangeVector].size();
	double numberComponents_mrs2 = mrs2.nodeRangeMaps[resolution][rangeVector].size();
	return min(numberComponents_mrs1, numberComponents_mrs2)/max(numberComponents_mrs1, numberComponents_mrs2);
}

//Compute the similarity between two nodes based on the their degrees (See Sections 4.3 and 5.4 in [1])
double vtkSimilarityBetweenMultiResolutionReebSpaces::SimilarityBetweenDegrees(int nodeId1, int nodeId2, int resolution){
	double degree1 = mrs1.jcnGraphs[resolution]->GetDegree(nodeId1);
	double degree2 = mrs2.jcnGraphs[resolution]->GetDegree(nodeId2);
	if(degree1 == 0 && degree2 == 0){
		return 1;
	}
	else{
		return (double)min(degree1, degree2)/(double)max(degree1, degree2);
	}
}

//Compute the similarity between two nodes (See Equation (3) in [1])
double vtkSimilarityBetweenMultiResolutionReebSpaces::SimilarityBetweenNodes(int nodeId1, int nodeId2, int resolution){
	return (weightVolume * SimilarityBetweenVolumes(nodeId1, nodeId2, resolution)) + (weightRangeMeasure * SimilarityBetweenRanges(nodeId1, nodeId2, resolution)) +(weightBettiZero * SimilarityBetweenNumberOfComponents(nodeId1, nodeId2, resolution)) +(weightDegree * SimilarityBetweenDegrees(nodeId1, nodeId2, resolution));
}

//Compute the similarity between the MPAIRS of two multi-resolution Reeb Spaces
double vtkSimilarityBetweenMultiResolutionReebSpaces::SimilarityBetweenMPAIRS(){
	long double finalSimilarity = 0.0;
	for (int resolution = 0; resolution < numberOfResolutions; ++resolution){
		double similarity_dimension = 0;
		vtkIntArray *sizeArray1 = vtkIntArray::SafeDownCast(mrs1.jcnGraphs[resolution]->GetVertexData()->GetArray("size"));
		vtkIntArray *sizeArray2 = vtkIntArray::SafeDownCast(mrs2.jcnGraphs[resolution]->GetVertexData()->GetArray("size"));
		for(int j = 0; j < MPAIRS[resolution].size(); j++){
			double size1 = sizeArray1->GetComponent(MPAIRS[resolution][j].first, 0);
			double size2 = sizeArray2->GetComponent(MPAIRS[resolution][j].second, 0);
			double weight = ((size1/totalSizeMrs1) + (size2/totalSizeMrs2))/2.0;
			similarity_dimension += weight * SimilarityBetweenNodes(MPAIRS[resolution][j].first, MPAIRS[resolution][j].second, resolution);
		}
		finalSimilarity += similarity_dimension;		
	}
	finalSimilarity = finalSimilarity/numberOfResolutions;
	return finalSimilarity;
}

//Check if the parents of the nodes to be matched form an MPAIR (See Section 5.1 in [1])
bool vtkSimilarityBetweenMultiResolutionReebSpaces::CheckParentMatching(int parentID1, int parentID2){
	int parentResolution = MPAIRS.size()-1;
	for(int i =0; i<MPAIRS[parentResolution].size(); i++){
		if(MPAIRS[parentResolution][i].first == parentID1 && MPAIRS[parentResolution][i].second == parentID2){
			return true;
		}
	} 
	return false;
}

//Check whether the mlists of the nodes to be matched are the same (topological consistency) - See Section 5.1 in [1]
bool vtkSimilarityBetweenMultiResolutionReebSpaces::CheckMLIST(vector<int> mlist1, vector<int> mlist2){
	int mlist1Size = mlist1.size();
	int mlist2Size = mlist2.size();
	if(mlist1Size != mlist2Size){
		return false;
	}
	for (int i = 0; i < mlist1Size; ++i){
		if(mlist1[i] != mlist2[i]){	
			return false;
		}
	}
	return true;
}

//Check the matching criteria (See Section 5.1 in [1])
bool vtkSimilarityBetweenMultiResolutionReebSpaces::CheckMatchingCriteria(int nodeId1, int nodeId2, int resolution){

	//nodeId1 and nodeId2 should not belong to any other MPAIR 
	if(IsMatchedMrs1[nodeId1] == true || IsMatchedMrs2[nodeId2] == true){
		return false;
	}
	//property : same range (already checked in the function FindCandidateMatchingNodes())

	//property : same matching parent
	if(resolution != 0){
		int parentId1 =  mrs1.jcnGraphs[resolution]->GetVertexData()->GetArray("parent")->GetComponent(nodeId1,0);
		int parentId2 =  mrs2.jcnGraphs[resolution]->GetVertexData()->GetArray("parent")->GetComponent(nodeId2,0);
		if(!CheckParentMatching(parentId1, parentId2)){		
			return false;
		}
	}
	//property : same MLISTs (topological consistency)
	vector<pair<vtkGraph*, int> > mdrgNodes1 = mdrgNodeMap1[nodeId1];
	vector<pair<vtkGraph*, int> > mdrgNodes2 = mdrgNodeMap2[nodeId2];
	for(int f = 0; f < this->Fields.size(); f++){
		vector<int> mlist1 = mlistMap1[mdrgNodes1[f]];
		vector<int> mlist2 = mlistMap2[mdrgNodes2[f]];
		if(!CheckMLIST(mlist1, mlist2)){
			return false;
		}
	}
	return true;
}

//Find candidate nodes satisfying the matching criteria (See Section 5.1 in [1])
vector<int> vtkSimilarityBetweenMultiResolutionReebSpaces::FindCandidateMatchingNodes(int nodeId, int resolution){
	vector<double> rangeValues;
	vtkDataSetAttributes *graphProperties = mrs1.jcnGraphs[resolution]->GetVertexData();
	for(int f = 0; f < this->Fields.size(); f++){
		rangeValues.push_back(graphProperties->GetArray(this->Fields[f])->GetComponent(nodeId,0));
	}
	vector <int>& nodeList = mrs2.nodeRangeMaps[resolution][rangeValues];//nodes in mrs2 which have the same range as nodeId in mrs1
	vector <int> candidateNodes;
	for (int i = 0; i < nodeList.size(); ++i){
		if(CheckMatchingCriteria(nodeId, nodeList[i],resolution)){
			candidateNodes.push_back(nodeList[i]);
		}
	}	
	return candidateNodes;
}

//Propogate the label along a branch in the monotonically increasing or decreasing direction (See Section 5.2 in [1])
void vtkSimilarityBetweenMultiResolutionReebSpaces::PropagateLabel(pair<vtkGraph*, int> node, int matchingLabel, bool isDecreasing, int reebSpace){
	vector<int> mlist;
	if(reebSpace == 1){
		mlist = mlistMap1[node];
		mlist.push_back(matchingLabel);
		mlistMap1[node] = mlist;
	}
	else{
		mlist = mlistMap2[node];
		mlist.push_back(matchingLabel);
		mlistMap2[node] = mlist;
	}
	vtkGraph* reebGraph = node.first;
	int nodeId = node.second;
	if(isDecreasing){
		int inDegree = reebGraph->GetInDegree(nodeId);//Number of adjacent nodes with lesser range
		if(inDegree == 1){
			vtkSmartPointer<vtkInEdgeIterator> edges = vtkSmartPointer<vtkInEdgeIterator>::New();
			reebGraph->GetInEdges(nodeId, edges);
			vtkInEdgeType inEdge= edges->Next();
			int nextNodeId = inEdge.Source;
			PropagateLabel(make_pair(reebGraph, nextNodeId), matchingLabel, isDecreasing, reebSpace);
		}
	}
	else{
		int outDegree = reebGraph->GetOutDegree(nodeId);//Number of adjacent nodes with greater range
		if(outDegree == 1){
			vtkSmartPointer<vtkOutEdgeIterator> edges = vtkSmartPointer<vtkOutEdgeIterator>::New();
			reebGraph->GetOutEdges(nodeId, edges);
			vtkOutEdgeType outEdge= edges->Next();
			int nextNodeId = outEdge.Target;
			PropagateLabel(make_pair(reebGraph, nextNodeId), matchingLabel, isDecreasing, reebSpace);
		}
	}
}

//Finding the MPAIRs (matching pairs of nodes) between JCNs of a particular resolution (see Algorithm 2 in [1])
void vtkSimilarityBetweenMultiResolutionReebSpaces::FindMatchingPairs(int resolution){
	vector< pair< int, int> > mpairsResolution;//Stores the pair of matched nodes for the JCNs of a particular resolution
	vector<int> matchingLabels(this->Fields.size(), -1);

	while(!queueMrs1Nodes.empty()){
		pair<double,int> highestSizeNode = queueMrs1Nodes.top();
		int nodeId = -highestSizeNode.second;
		queueMrs1Nodes.pop();
		vector <int> candidateNodes = FindCandidateMatchingNodes(nodeId, resolution);
		if(candidateNodes.size() > 0){
			//finding a single matching node from candidate nodes which has maximum matching with the chosen node
			double matchingValue = SimilarityBetweenNodes(nodeId, candidateNodes[0], resolution);
			int matchingNodeId=candidateNodes[0];
			for (int i = 1; i < candidateNodes.size(); ++i){
				double currentMatchingValue = SimilarityBetweenNodes(nodeId, candidateNodes[i], resolution);
				if(currentMatchingValue > matchingValue ){
					matchingValue = currentMatchingValue;
					matchingNodeId = candidateNodes[i];
				}
			}
			//create the MPAIR
			mpairsResolution.push_back(make_pair(nodeId,matchingNodeId));
			IsMatchedMrs1[nodeId] = true;
			IsMatchedMrs2[matchingNodeId] = true;
			vector<pair<vtkGraph*, int> >	mdrgNodeList1 = mdrgNodeMap1[nodeId];
			vector<pair<vtkGraph*, int> >	mdrgNodeList2 = mdrgNodeMap2[matchingNodeId];
			for(int f = 0; f < this->Fields.size(); f++){
				//update matching label for field f and propagate the label
				matchingLabels[f] = matchingLabels[f] + 1;
				pair<vtkGraph*, int> node = mdrgNodeList1[f];
				vtkGraph* reebGraph = node.first;
				int mdrgNodeId = node.second;
				PropagateLabel(node, matchingLabels[f], true, 1);
				int outDegree = reebGraph->GetOutDegree(mdrgNodeId);//Number of adjacent nodes with greater range
				if(outDegree == 1){
					vtkSmartPointer<vtkOutEdgeIterator> edges = vtkSmartPointer<vtkOutEdgeIterator>::New();
					reebGraph->GetOutEdges(mdrgNodeId, edges);
					vtkOutEdgeType outEdge = edges->Next();
					int nextNodeId = outEdge.Target;
					PropagateLabel(make_pair(reebGraph, nextNodeId), matchingLabels[f], false, 1);
				}

				node = mdrgNodeList2[f];
				reebGraph = node.first;
				mdrgNodeId = node.second;
				PropagateLabel(node, matchingLabels[f], true, 2);
				outDegree = reebGraph->GetOutDegree(mdrgNodeId);//Number of adjacent nodes with greater range
				if(outDegree == 1){
					vtkSmartPointer<vtkOutEdgeIterator> edges = vtkSmartPointer<vtkOutEdgeIterator>::New();
					reebGraph->GetOutEdges(mdrgNodeId, edges);
					vtkOutEdgeType outEdge = edges->Next();
					int nextNodeId = outEdge.Target;
					PropagateLabel(make_pair(reebGraph, nextNodeId), matchingLabels[f], false, 2);
				}
			}
		}
	}
	//Updating MPAIRs 
	this->MPAIRS.push_back(mpairsResolution);
}

//Add the parameters for computing the two multi-resolution Reeb spaces and the similarity between them
void vtkSimilarityBetweenMultiResolutionReebSpaces::AddParameters(char *fieldName, double minimumValue,double maximumValue,double slabwidth){
	char *f = new char[strlen(fieldName) + 1];
	strcpy(f, fieldName);
	this->Fields.push_back(f);
	this->minimumValues.push_back(minimumValue);
	this->maximumValues.push_back(maximumValue);
	this->slabwidths.push_back(slabwidth);
}

//Compute the total size of the fiber-components corresponding to all the nodes in each of the multi-resolution Reeb spaces
void vtkSimilarityBetweenMultiResolutionReebSpaces::ComputeTotalSize(){
	//Compute the total size of the nodes in MRS 1
	this->totalSizeMrs1 = 0;
	int nrVertices = mrs1.jcnGraphs[0]->GetNumberOfVertices();
	vtkIntArray *sizeArray = vtkIntArray::SafeDownCast(mrs1.jcnGraphs[0]->GetVertexData()->GetArray("size"));
	for(int i = 0; i < nrVertices; i++){
		totalSizeMrs1 = totalSizeMrs1 + sizeArray->GetComponent(i,0);
	}
	//Compute the total size of the nodes in MRS 2
	this->totalSizeMrs2 = 0;
	nrVertices = mrs2.jcnGraphs[0]->GetNumberOfVertices();
	sizeArray = vtkIntArray::SafeDownCast(mrs2.jcnGraphs[0]->GetVertexData()->GetArray("size"));
	for(int i = 0; i < nrVertices; i++){
		totalSizeMrs2 = totalSizeMrs2 + sizeArray->GetComponent(i,0);
	}
}

//Compute the similarity between two multi-resolution Reeb spaces with finest resolution JCNs "JCN1" and "JCN2", respectively.
double vtkSimilarityBetweenMultiResolutionReebSpaces::ComputeSimilarity(vtkGraph * JCN1, vtkGraph * JCN2, int numberOfResolutions){
	this->numberOfResolutions = numberOfResolutions;
	mrs1.ConstructMrs(JCN1, Fields, minimumValues, maximumValues, slabwidths, numberOfResolutions);
	mrs2.ConstructMrs(JCN2, Fields, minimumValues, maximumValues, slabwidths, numberOfResolutions);
	ComputeTotalSize();
	int numberOfNodesInMRS1=0;
	int numberOfNodesInMRS2=0;
	int numberOfMPAIRS=0;
	int volumeOfMatchedNodesInMRS1 = 0,volumeOfMatchedNodesInMRS2 = 0;
	for(int resolution =0; resolution < this->numberOfResolutions; resolution++){
		Initialize(resolution);
		FindMatchingPairs(resolution);//finding the MPAIRs between JCNs of a particular resolution
		numberOfNodesInMRS1 = numberOfNodesInMRS1 + mrs1.jcnGraphs[resolution]->GetNumberOfVertices();
		numberOfNodesInMRS2 = numberOfNodesInMRS2 + mrs2.jcnGraphs[resolution]->GetNumberOfVertices();
		numberOfMPAIRS = numberOfMPAIRS + this->MPAIRS[resolution].size();
		for(int i = 0; i <  this->MPAIRS[resolution].size(); i++){
			volumeOfMatchedNodesInMRS1 = volumeOfMatchedNodesInMRS1 + mrs1.jcnGraphs[resolution]->GetVertexData()->GetArray("size")->GetComponent(this->MPAIRS[resolution][i].first,0);
			volumeOfMatchedNodesInMRS2 = volumeOfMatchedNodesInMRS2 + mrs2.jcnGraphs[resolution]->GetVertexData()->GetArray("size")->GetComponent(this->MPAIRS[resolution][i].second,0);
		}
	}
	return SimilarityBetweenMPAIRS();
}
