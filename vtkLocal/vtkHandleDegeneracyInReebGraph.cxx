/*=========================================================================

 *	File: vtkHandleDegeneracyInReebGraph.cxx
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
//Reference:
//[1] J. Tu, M. Hajij, and P. Rosen, “Propagate and pair: A single-pass approach to critical point pairing in reeb graphs,” in Advances in Visual Computing,
// Eds. Springer International Publishing, 2019.

#include "vtkHandleDegeneracyInReebGraph.h"
#include "vtkAdjacentVertexIterator.h"
#include "vtkObjectFactory.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkInEdgeIterator.h"
#include "vtkOutEdgeIterator.h"
#include<bits/stdc++.h>
using namespace std;

vtkStandardNewMacro(vtkHandleDegeneracyInReebGraph);


#define VTK_CREATE(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

vtkHandleDegeneracyInReebGraph::vtkHandleDegeneracyInReebGraph(){

}

vtkHandleDegeneracyInReebGraph::~vtkHandleDegeneracyInReebGraph(){
  
}

//Sets the Reeb graph and the name of the associate scalar-field
void vtkHandleDegeneracyInReebGraph::SetParameters(vtkReebGraph* reebGraph, char *fieldName){
	char *f = new char[strlen(fieldName) + 1];
	strcpy(f, fieldName);
	this->fieldName = f;
	this->reebGraph = reebGraph;
}

//Breaks a non-degenerate minimum into a minimum and an up-fork
void vtkHandleDegeneracyInReebGraph::HandleDegenerateMinimum(int nodeId){
	//Creating a duplicate vertex in the reeb graph
	int duplicateVertexId = reebGraph->GetNumberOfVertices();
	reebGraph->AddVertex();
	reebGraph->AddEdge(duplicateVertexId, nodeId);//Connecting the vertex with its duplicate
	vtkSmartPointer<vtkDoubleArray> field;
	field = vtkDoubleArray::SafeDownCast(reebGraph->GetVertexData()->GetArray(fieldName));
	if((reebGraph->GetNumberOfVertices() - 1) < field->GetNumberOfTuples()){
		field->SetComponent(reebGraph->GetNumberOfVertices()-1, 0, field->GetComponent(nodeId,0) - epsilon);
	}
	else{
		field->InsertNextTuple1(field->GetComponent(nodeId,0) - epsilon);//Adding the new field value to the duplicate vertex
	}
}

//Breaks a non-degenerate maximum into a maximum and a down-fork
void vtkHandleDegeneracyInReebGraph::HandleDegenerateMaximum(int nodeId){
	//Creating a duplicate vertex in the reeb graph
	int duplicateVertexId = reebGraph->GetNumberOfVertices();
	reebGraph->AddVertex();
	reebGraph->AddEdge(nodeId, duplicateVertexId);//Connecting the vertex with its duplicate
	vtkSmartPointer<vtkDoubleArray> field;
	field = vtkDoubleArray::SafeDownCast(reebGraph->GetVertexData()->GetArray(fieldName));
	if((reebGraph->GetNumberOfVertices() - 1) < field->GetNumberOfTuples()){
		field->SetComponent(reebGraph->GetNumberOfVertices()-1, 0, field->GetComponent(nodeId,0) + epsilon);
	}
	else{
		field->InsertNextTuple1(field->GetComponent(nodeId,0) + epsilon);//Adding the new field value to the duplicate vertex
	}
}

//Break a degenerate up-fork with up-degree n into (n-1) up-forks
void vtkHandleDegeneracyInReebGraph::HandleDegenerateUpFork(int nodeId){
	int nrVertices = reebGraph->GetNumberOfVertices();
	int outDegree = reebGraph->GetOutDegree(nodeId);//Number of adjacent nodes with greater range value
	vtkSmartPointer<vtkOutEdgeIterator> edges = vtkSmartPointer<vtkOutEdgeIterator>::New();
	reebGraph->GetOutEdges(nodeId, edges);
	vector<int> adjacentVertices;
	while(edges->HasNext()){
		adjacentVertices.push_back(edges->Next().Target);
	}
	for(int i = 0; i < (outDegree - 2);i++){
		int newDuplicateVertexId = reebGraph->GetNumberOfVertices();
		reebGraph->AddVertex();
		vtkSmartPointer<vtkDoubleArray> field;
		field = vtkDoubleArray::SafeDownCast(reebGraph->GetVertexData()->GetArray(fieldName));
		if((reebGraph->GetNumberOfVertices() - 1) < field->GetNumberOfTuples()){
			field->SetComponent(reebGraph->GetNumberOfVertices()-1, 0, field->GetComponent(nodeId,0)  + ((i+1)*epsilon));
		}
		else{
			field->InsertNextTuple1(field->GetComponent(nodeId,0) + ((i+1)*epsilon));//Adding the new field value to the duplicate vertex
		}

		if(i==0){
			reebGraph->AddEdge(nodeId, newDuplicateVertexId);//Connecting the vertex with its duplicate
		}
		else{
			reebGraph->AddEdge(newDuplicateVertexId - 1, newDuplicateVertexId);//Connecting the vertex with its duplicate			
		}
		int targetNodeId = adjacentVertices[i];
		reebGraph->RemoveEdge(reebGraph->GetEdgeId(nodeId, adjacentVertices[i]));
		reebGraph->AddEdge(newDuplicateVertexId, targetNodeId);//Changing the edge from nodeId to an edge from duplicateNodeId
		if(i==(outDegree - 3)){//The upper-most duplicate vertex has two upward edges. Adding the second upward edge.
			reebGraph->RemoveEdge(reebGraph->GetEdgeId(nodeId, adjacentVertices[i+1]));
			reebGraph->AddEdge(newDuplicateVertexId, adjacentVertices[i+1]);//Changing the edge from nodeId to an edge from duplicateNodeId		
		}
	}
}

//Break a degenerate down-fork with up-degree n into (n-1) down-forks
void vtkHandleDegeneracyInReebGraph::HandleDegenerateDownFork(int nodeId){
	int nrVertices = reebGraph->GetNumberOfVertices();
	int inDegree = reebGraph->GetInDegree(nodeId);//Number of adjacent nodes with lesser range value
	vtkSmartPointer<vtkInEdgeIterator> edges = vtkSmartPointer<vtkInEdgeIterator>::New();
	reebGraph->GetInEdges(nodeId, edges);
	vector<int> adjacentVertices;
	while(edges->HasNext()){
		adjacentVertices.push_back(edges->Next().Source);
	}
	for(int i = 0; i < (inDegree - 2); i++){
		int newDuplicateVertexId = reebGraph->GetNumberOfVertices();
		reebGraph->AddVertex();
		vtkSmartPointer<vtkDoubleArray> field;
		field = vtkDoubleArray::SafeDownCast(reebGraph->GetVertexData()->GetArray(fieldName));
		if((reebGraph->GetNumberOfVertices() - 1) < field->GetNumberOfTuples()){
			field->SetComponent(reebGraph->GetNumberOfVertices()-1, 0, field->GetComponent(nodeId,0)  - ((i+1)*epsilon));
		}
		else{
			field->InsertNextTuple1(field->GetComponent(nodeId,0) - ((i+1)*epsilon));//Adding the new field value to the duplicate vertex
		}
		if(i==0){
			reebGraph->AddEdge(newDuplicateVertexId,nodeId);//Connecting the vertex with its duplicate
		}
		else{
			reebGraph->AddEdge(newDuplicateVertexId,newDuplicateVertexId - 1);//Connecting the vertex with its duplicate			
		}
		reebGraph->RemoveEdge(reebGraph->GetEdgeId(nodeId, adjacentVertices[i]));
		reebGraph->AddEdge(adjacentVertices[i],newDuplicateVertexId);//Changing the edge from nodeId to an edge from duplicateNodeId
		if(i == (inDegree - 3)){//The lower-most duplicate vertex has two downward edges. Adding the second downward edge.
			reebGraph->RemoveEdge(reebGraph->GetEdgeId(adjacentVertices[i+1], nodeId));
			reebGraph->AddEdge(adjacentVertices[i+1], newDuplicateVertexId);//Changing the edge from nodeId to an edge from duplicateNodeId		
		}
	}
}

//Break a double-fork into a down-fork and an up-fork
void vtkHandleDegeneracyInReebGraph::HandleDegenerateDoubleFork(int nodeId){
	int newDuplicateVertexId = reebGraph->GetNumberOfVertices();
	reebGraph->AddVertex();
	vtkSmartPointer<vtkDoubleArray> field;
	field = vtkDoubleArray::SafeDownCast(reebGraph->GetVertexData()->GetArray(fieldName));
	if((reebGraph->GetNumberOfVertices() - 1) < field->GetNumberOfTuples()){
		field->SetComponent(reebGraph->GetNumberOfVertices()-1, 0, field->GetComponent(nodeId,0)  - epsilon);
	}
	else{
		field->InsertNextTuple1(field->GetComponent(nodeId,0) - epsilon);//Adding the new field value to the duplicate vertex
	}


	vtkSmartPointer<vtkInEdgeIterator> edges = vtkSmartPointer<vtkInEdgeIterator>::New();
	reebGraph->GetInEdges(nodeId, edges);
	vector<int> adjacentVertices;
	while(edges->HasNext()){
		adjacentVertices.push_back(edges->Next().Source);
	}
	for(int i = 0; i < adjacentVertices.size(); i++){
		reebGraph->RemoveEdge(reebGraph->GetEdgeId(adjacentVertices[i], nodeId));
		reebGraph->AddEdge(adjacentVertices[i],newDuplicateVertexId);//Changing the edge from nodeId to an edge from duplicateNodeId
	}
	reebGraph->AddEdge(newDuplicateVertexId, nodeId);//Connecting the vertex with its duplicate
}

/*Removes degeneracy in the Reeb graph by:
1.  Handling degenerate critical nodes of the Reeb graph which violate the first Morse condition
2.  Ensuring that no two critical nodes of the Reeb graph have the same critical value (second Morse condition)
*/
void vtkHandleDegeneracyInReebGraph::HandleDegeneracy(){
	int nrVertices = reebGraph->GetNumberOfVertices();
	int inDegree,outDegree;
	vtkSmartPointer<vtkDataArray> field;
	field = reebGraph->GetVertexData()->GetArray(fieldName);
	int nrComponents = field->GetNumberOfTuples();

	for(int j = 0; j < nrVertices; j++){
		inDegree = reebGraph->GetInDegree(j);
		outDegree = reebGraph->GetOutDegree(j);
		if(inDegree == 0 and outDegree >= 2){
			HandleDegenerateMinimum(j);
		}
		if(outDegree == 0 and inDegree >= 2){
			HandleDegenerateMaximum(j);
		}
	}
	for(int j = 0; j < nrVertices; j++){
		inDegree = reebGraph->GetInDegree(j);
		outDegree = reebGraph->GetOutDegree(j);
		if(inDegree >= 2 and outDegree >= 2){
			HandleDegenerateDoubleFork(j);
		}
	}
	for(int j = 0; j < nrVertices; j++){
		inDegree = reebGraph->GetInDegree(j);
		outDegree = reebGraph->GetOutDegree(j);
		if(outDegree >= 3){
			HandleDegenerateUpFork(j);
		}
		else if(inDegree >= 3){
			HandleDegenerateDownFork(j);
		}
	}
	AssignUniqueRangeValuesToCriticalNodes();
}

//Ensures that no two critical nodes of the Reeb graph have the same critical value (second Morse condition)
void vtkHandleDegeneracyInReebGraph::AssignUniqueRangeValuesToCriticalNodes(){
	set<double> rangeValues;//Stores the range values in non-decreasing order
	map<double, vector<int>> criticalRangeNodesMap;//Store the node indices corresponding to a critical value
	int nrVertices = reebGraph->GetNumberOfVertices();
	vtkSmartPointer<vtkDoubleArray> field;
	field = vtkDoubleArray::SafeDownCast(reebGraph->GetVertexData()->GetArray(fieldName));
    	for(int j = 0; j < nrVertices; j++){
		double range = field->GetComponent(j, 0);
    		rangeValues.insert(field->GetComponent(j, 0));
    		int inDegree = reebGraph->GetInDegree(j);
		int outDegree = reebGraph->GetOutDegree(j);
		if(inDegree != 1 || outDegree != 1){
			if(criticalRangeNodesMap.find(range) == criticalRangeNodesMap.end()){
				vector<int> v;
				v.push_back(j);
				criticalRangeNodesMap[range] = v;
			}
			else{
				criticalRangeNodesMap[range].push_back(j);
			}
		}
    	}
    	for(std::map<double, vector<int>>::iterator iter = criticalRangeNodesMap.begin(); iter != criticalRangeNodesMap.end(); ++iter){
		double range =  iter->first;
		vector<int> criticalNodes = iter->second;
		if(criticalNodes.size() > 1){
			double nextRange;
			int numberOfCriticalNodes = criticalNodes.size();
			set<double>::iterator itr;
			for (itr = rangeValues.begin(); *itr != range; itr++){
			}
			itr++;
			if(itr != rangeValues.end()){
				nextRange = *itr;
			}
			else{
				nextRange = range + epsilon;
			}
			double inc = (nextRange - range)/double(max(default_number_of_range_values,numberOfCriticalNodes));
			for(int i = 1; i < criticalNodes.size(); i++){
				field->SetComponent(criticalNodes[i], 0, range + (inc * i));//Assigning an unique range value to each critical node	
			}
		}
	}
}
