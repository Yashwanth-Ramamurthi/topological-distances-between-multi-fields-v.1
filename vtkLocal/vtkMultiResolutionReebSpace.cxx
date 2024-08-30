/*=========================================================================

 *	File: vtkMultiResolutionReebSpace.cxx
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
//Reference:
//[1] Yashwanth Ramamurthi and Amit Chattopadhyay. "A Topological Similarity Measure Between Multi-Resolution Reeb Spaces," in IEEE Transactions on Visualization and Computer Graphics, vol. 28, no. 12, pp. 4360-4374, 1 Dec. 2022

#include "vtkMultiResolutionReebSpace.h"
#include "vtkMultiDimensionalReebGraph.h"
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
#include <math.h>
#include <fstream>
#include <vtkVariant.h>
#include <vtkVariantArray.h>
using namespace std;

vtkStandardNewMacro(vtkMultiResolutionReebSpace);

#define VTK_CREATE(type,name) \
vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

vtkMultiResolutionReebSpace::vtkMultiResolutionReebSpace(){
}

vtkMultiResolutionReebSpace::~vtkMultiResolutionReebSpace(){
}

void vtkMultiResolutionReebSpace::PrintSelf(ostream& os, vtkIndent indent){
	this->Superclass::PrintSelf(os, indent);
}

//Find the nodes in the finer resolution JCN which will be merged to a single node in the coarser resolution JCN
//Nodes in "component" (returned by this procedure) will be merged to a single node in the coarser resolution JCN
vector<int> vtkMultiResolutionReebSpace::UnionFind(vtkGraph *jcn, int nodeId, vector<bool> &visited,vector<int> component, vector<vector<pair<double,double> > > rangePairs){
	visited[nodeId] = true;
	component.push_back(nodeId);
	vtkSmartPointer<vtkAdjacentVertexIterator> adj = vtkSmartPointer<vtkAdjacentVertexIterator>::New();
	jcn->GetAdjacentVertices(nodeId, adj);
	while (adj->HasNext()){
		int adjacentNodeId = adj->Next();
		if(visited[adjacentNodeId] == false){
			bool check = true;
			for(int f = 0; f<this->Fields.size() && check==true; f++){//check will be true at the end of this for loop if the nodes nodeId and adjacentNodeId should be merged
				check = false;
				double sourceNodeRange = jcn->GetVertexData()->GetArray(this->Fields[f])->GetComponent(nodeId, 0);
				double targetNodeRange = jcn->GetVertexData()->GetArray(this->Fields[f])->GetComponent(adjacentNodeId, 0);
				for(int i = 0; i < rangePairs[f].size(); i++){
					pair<double,double> rangePair = rangePairs[f][i];
					if((rangePair.first==sourceNodeRange || rangePair.second==sourceNodeRange) && (rangePair.first==targetNodeRange || rangePair.second==targetNodeRange)){
						check = true;
						break;
					}
				}
			}
			if(check == true){
				component = UnionFind(jcn, adjacentNodeId, visited, component, rangePairs);
			}
		}
	}
	return component;
}


//Obtain the id of the component containing the node with id "nodeId"
int vtkMultiResolutionReebSpace::GetComponentNumber(int nodeId, vector<vector<int> > components){
	for(int i = 0; i < components.size(); i++){
		for(int j = 0; j < components[i].size(); j++){
			if(components[i][j] == nodeId){
				return i;
			}
		}
	}
}

//Create a coarser JCN from a fine JCN. Nodes in the finer JCN will be merged to obtain the coarser JCN
vtkSmartPointer<vtkUndirectedGraph> vtkMultiResolutionReebSpace::ComputeCoarseJCN(vtkGraph *fineJCN,vector<double> minimumValues, vector<double> maximumValues, vector<double> slabwidth){
	int nrEdges = fineJCN->GetNumberOfEdges();
	int nrVertices = fineJCN->GetNumberOfVertices();
	vector<vector<pair<double,double> > > rangePairs;
	// creating pairs of ranges that needs to be merged together to create coarser quantized ranges
	for(int f = 0; f < this->Fields.size(); f++){
		vector<pair<double,double> > rangePairsForFieldF; //Stores the pairs of ranges for field f in the present resolution which will merge into a single range in the coarser resolution
		for (double value = minimumValues[f]; value < maximumValues[f]; value += 2*slabwidth[f]){
			pair<double,double> pair_set;
			pair_set.first = value;
			pair_set.second = value + slabwidth[f];
			rangePairsForFieldF.push_back(pair_set);
		}
		rangePairs.push_back(rangePairsForFieldF);
	}
	//Performing UnionFind to find out the nodes which will be merged
	vector<bool> visited(nrVertices,false);
	 vector<vector<int> > components;
	for(int i = 0; i < nrVertices; i++){
		if(visited[i] == false){
			vector<int> component;
			components.push_back( UnionFind(fineJCN, i, visited, component, rangePairs) );	//Nodes in "component" will be merged to a single node in the coarse JCN
		}		
	}
	//creating the coarse JCN with the with information about vertexData, edgeData, size, field values.
	int nrVerticesCoarseJCN=components.size();
	VTK_CREATE(vtkMutableGraphHelper, builder);
	VTK_CREATE(vtkMutableUndirectedGraph, coarseJCN);
	builder->SetGraph(coarseJCN);
	//Adding vertex Id
	for (int i=0; i<nrVerticesCoarseJCN; i++){
		builder->AddVertex();
	}

	vector<pair<int, int> > edges;
	for(int i=0;i<nrEdges;i++){
		vtkIdType sourceOfEdge = fineJCN->GetSourceVertex(i);
		vtkIdType targetOfEdge = fineJCN->GetTargetVertex(i);
		int sourceCompId = GetComponentNumber(sourceOfEdge, components);
		int targetCompId = GetComponentNumber(targetOfEdge, components);
		if(sourceCompId != targetCompId){
			pair<int,int> edge = make_pair(sourceCompId, targetCompId);
			if(find(edges.begin(),edges.end(),make_pair(sourceCompId, targetCompId)) == edges.end() && find(edges.begin(),edges.end(),make_pair(targetCompId, sourceCompId)) == edges.end()){
				builder->AddEdge(sourceCompId,targetCompId);
				edges.push_back(edge);
			}
		}
	}
	//Store the size for each node in the coarse JCN
	VTK_CREATE(vtkIntArray, nodeSize);
	nodeSize->SetName("size");
	nodeSize->SetNumberOfValues(nrVerticesCoarseJCN);
	for (int i=0; i<nrVerticesCoarseJCN; i++){//Size of a node X in the coarser JCN is the sum of the sizes of the nodes in the finer JCN which are merged to form X
		int size = 0;
		for(int c = 0; c < components[i].size(); c++){
			size = size + fineJCN->GetVertexData()->GetArray("size")->GetComponent(components[i][c],0);
		}
		nodeSize->SetValue(i, size);
	}
	coarseJCN->GetVertexData()->AddArray(nodeSize);
	//Store the field values for each node in the coarseJCN
	for(int f = 0; f < this->Fields.size(); f++){
		VTK_CREATE(vtkFloatArray, field);
		field->SetName(this->Fields[f]);
		field->SetNumberOfValues(nrVerticesCoarseJCN);
		for(int i = 0; i < nrVerticesCoarseJCN; i++){
			int nodeInFineJCN = components[i][0];
			double range = fineJCN->GetVertexData()->GetArray(this->Fields[f])->GetComponent(nodeInFineJCN, 0);
			if((int)((range - minimumValues[f])/slabwidth[f])%2 == 0){
				field->SetValue(i, range);
			}
			else{
				field->SetValue(i, range - slabwidth[f]);
			}
		}
		coarseJCN->GetVertexData()->AddArray(field);		
	}

	//Store parent-child relationship between the nodes of the coarser and finer JCNs
	vtkSmartPointer<vtkVariantArray> children = vtkSmartPointer<vtkVariantArray>::New();
	children->SetName("children");
	VTK_CREATE(vtkIdTypeArray, parent);
	parent->SetName("parent");
	parent->SetNumberOfValues(nrVertices);
	for(int i = 0; i < nrVerticesCoarseJCN; i++){
		VTK_CREATE(vtkIdTypeArray, childrenOfANode);
		childrenOfANode->SetName("childrenOfANode");
		childrenOfANode->SetNumberOfValues(components[i].size());
		for(int j = 0; j < components[i].size(); j++){
			childrenOfANode->SetValue(j,components[i][j]);
			parent->SetValue(components[i][j],i);
		}
		children->InsertNextValue(vtkVariant(childrenOfANode));
	}
	coarseJCN->GetVertexData()->AddArray(children);
	fineJCN->GetVertexData()->AddArray(parent);
	VTK_CREATE(vtkIdTypeArray, nodeIds);
	nodeIds->SetName("jcnids");
	nodeIds->SetNumberOfValues(nrVerticesCoarseJCN);
	for (vtkIdType i = 0; i < nrVerticesCoarseJCN; i++)
	{
	     nodeIds->SetValue(i, i);//copying the origical JCN-ids  
    	}
 	coarseJCN->GetVertexData()->AddArray(nodeIds);//Adding JCN-ids in Input Graph
	
	return coarseJCN;
}

//Construct the multi-resolution Reeb Space by taking the JCN in the finest resolution as input
void vtkMultiResolutionReebSpace::ConstructMrs(vtkGraph* finestResolutionJCN, vector <char*> fields, vector<double> minimumValues, vector<double> maximumValues, vector<double> slabwidth, int numberOfResolutions){
	this->numberOfResolutions = numberOfResolutions;
	for(int f = 0; f < fields.size(); f++){
		this->Fields.push_back(fields[f]);
	}
	vtkSmartPointer<vtkUndirectedGraph> fineJCN = vtkUndirectedGraph::SafeDownCast(finestResolutionJCN);
	int nrVertices = fineJCN->GetNumberOfVertices();
	VTK_CREATE(vtkIdTypeArray, nodeIds);
	nodeIds->SetName("jcnids");
	nodeIds->SetNumberOfValues(nrVertices);
	for (vtkIdType i = 0; i < nrVertices; i++)
	{
	     nodeIds->SetValue(i, i);//copying the origical JCN-ids  
    	}
 	fineJCN->GetVertexData()->AddArray(nodeIds);//Adding JCN-ids in Input Graph
	vtkSmartPointer<vtkUndirectedGraph> coarseJCN;
	jcnGraphs.push_back(fineJCN);
	for(int resolution = 1; resolution  < numberOfResolutions; resolution++){//Create coarser JCNs, Nodes in fineJCN will be merged to obtain coarseJCN
		coarseJCN = ComputeCoarseJCN(fineJCN, minimumValues, maximumValues, slabwidth);
		jcnGraphs.insert(jcnGraphs.begin(), coarseJCN);
		fineJCN=coarseJCN;

		//The size of a slab in the coarse JCN is two times the size in the fine JCN
		for(int f = 0; f < this->Fields.size(); f++){
			slabwidth[f]*=2.0;
		}
	}
	for(int r = 0; r < numberOfResolutions; r++){
		vtkGraph *jcn = jcnGraphs[r];
		map<vector<double>, vector<int> >  nodeRangeMapCurrentResolution;//Stores the list of nodes corresponding to every combination of range values of the fields
		vector<double> minimumRangesCurrentResolution;//Stores the minimum value for each field
		vector<double> maximumRangesCurrentResolution;//Stores the maximum value for each field
		for(int f = 0; f < this->Fields.size(); f++){
			minimumRangesCurrentResolution.push_back(jcn->GetVertexData()->GetArray(this->Fields[f])->GetComponent(0, 0));
			maximumRangesCurrentResolution.push_back(jcn->GetVertexData()->GetArray(this->Fields[f])->GetComponent(0, 0));
		}
		for(int i = 0; i < jcn->GetNumberOfVertices(); i++){
			vector<double> rangeVector;
			for(int f = 0; f < this->Fields.size(); f++){
				double range = jcn->GetVertexData()->GetArray(this->Fields[f])->GetComponent(i, 0);
				rangeVector.push_back(range);
				if(range < minimumRangesCurrentResolution[f])
					minimumRangesCurrentResolution[f] = range;
				if(range > maximumRangesCurrentResolution[f])
					maximumRangesCurrentResolution[f] = range;
			}
			if(nodeRangeMapCurrentResolution.find(rangeVector) == nodeRangeMapCurrentResolution.end()){
				vector<int> nodeList;
				nodeList.push_back(i);
				nodeRangeMapCurrentResolution[rangeVector] = nodeList;
			}
			else{
				vector<int> nodeList = nodeRangeMapCurrentResolution[rangeVector];
				nodeList.push_back(i);
				nodeRangeMapCurrentResolution[rangeVector] = nodeList;
			}	
		}
		nodeRangeMaps.push_back(nodeRangeMapCurrentResolution);
		minimumRanges.push_back(minimumRangesCurrentResolution);
		maximumRanges.push_back(maximumRangesCurrentResolution);
	}
}
