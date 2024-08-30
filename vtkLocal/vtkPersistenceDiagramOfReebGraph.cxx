/*=========================================================================

 *	File: vtkPersistenceDiagramOfReebGraph.cxx
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
//Reference:
//Ulrich Bauer, Xiaoyin Ge, and Yusu Wang. 2014. Measuring Distance between Reeb Graphs. In Proceedings of the thirtieth annual symposium on Computational geometry (SOCG'14).
//Association for Computing Machinery, New York, NY, USA

//See also the following article for the construction of different persistence diagrams:
//Yashwanth Ramamurthi and Amit Chattopadhyay, "A Topological Distance Between Multi-Fields Based on Multi-Dimensional Persistence Diagrams," in IEEE Transactions on Visualization and Computer Graphics, vol. 30, no. 9, pp. 5939-5952, Sept. 2024

#include "vtkPersistenceDiagramOfReebGraph.h"
#include "vtkObjectFactory.h"
#include "vtkDoubleArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkInEdgeIterator.h"
#include "vtkOutEdgeIterator.h"
#include "vtkIdTypeArray.h"

vtkStandardNewMacro(vtkPersistenceDiagramOfReebGraph);

#define VTK_CREATE(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

vtkPersistenceDiagramOfReebGraph::vtkPersistenceDiagramOfReebGraph(){

}

vtkPersistenceDiagramOfReebGraph::~vtkPersistenceDiagramOfReebGraph(){
}

//Obtain the persistence diagram of a specific type
//type 0 => 0-dimensional Persistence Diagram corresponding to the sublevelset filtration 
//type 1 => 0-dimensional Persistence Diagram corresponding to the superlevelset filtration
//type 2 => 0-dimensional Extended Persistence Diagram
//type 3 => 1-dimensional Extended Persistence Diagram
vector<pair<double, double> > vtkPersistenceDiagramOfReebGraph::GetPersistenceDiagram(int type){
	if(type == 0){
		return persistenceDiagramSublevelSetFiltration;
	}
	else if(type == 1){
		return persistenceDiagramSuperLevelSetFiltration;
	}
	else if(type == 2){
		return extendedPersistenceDiagram_dimension_0;	
	}
	else if(type == 3){
		return extendedPersistenceDiagram_dimension_1;
	}
}

//Returns the id of the node with the minimum range (f) value by traversing the Reeb graph from nodeId to all the nodes with range (f) value  at most "maximumRangeValue"
vtkIdType vtkPersistenceDiagramOfReebGraph::GetMinimumNode(int nodeId, double maximumRangeValue, vector<int> &path){
	path.push_back(nodeId);
	if(reebGraph->GetInDegree(nodeId) == 0){
		return nodeId;
	}
	vector<int> adjacentVertices;
	vtkSmartPointer<vtkDoubleArray> field;
	field = vtkDoubleArray::SafeDownCast(reebGraph->GetVertexData()->GetArray(fieldName));
	vtkSmartPointer<vtkInEdgeIterator> inEdges = vtkSmartPointer<vtkInEdgeIterator>::New();
	reebGraph->GetInEdges(nodeId, inEdges);
	while(inEdges->HasNext()){
		adjacentVertices.push_back(inEdges->Next().Source);
	}
	vtkSmartPointer<vtkOutEdgeIterator> outEdges = vtkSmartPointer<vtkOutEdgeIterator>::New();
	reebGraph->GetOutEdges(nodeId, outEdges);
	while(outEdges->HasNext()){
		adjacentVertices.push_back(outEdges->Next().Target);
	}
	int minimumNode = -1;
	for(int i = 0; i < adjacentVertices.size(); i++){
		int adjacentnodeId = adjacentVertices[i];
		double adjacentNodeRange = field->GetComponent(adjacentnodeId, 0);
		if(find(path.begin(), path.end(), adjacentnodeId) == path.end() && adjacentNodeRange < maximumRangeValue){
			int currentMinimumNode = GetMinimumNode(adjacentnodeId, maximumRangeValue, path);
			if(minimumNode == -1){
				minimumNode = currentMinimumNode;
			}
			else if(currentMinimumNode != -1){
				if(field->GetComponent(currentMinimumNode, 0) < field->GetComponent(minimumNode, 0) ){
					minimumNode = currentMinimumNode;
				}
			}
		}
	}
	return minimumNode;
}

//Returns the id of the node with the maximum range (f) value by traversing the Reeb graph from nodeId to all the nodes with range (f) value  at least "minimumRangeValue"
vtkIdType vtkPersistenceDiagramOfReebGraph::GetMaximumNode(int nodeId, double minimumRangeValue, vector<int> &path){
	path.push_back(nodeId);
	if(reebGraph->GetOutDegree(nodeId) == 0){
		return nodeId;
	}
	vector<int> adjacentVertices;
	vtkSmartPointer<vtkDoubleArray> field;
	field = vtkDoubleArray::SafeDownCast(reebGraph->GetVertexData()->GetArray(fieldName));
	vtkSmartPointer<vtkInEdgeIterator> inEdges = vtkSmartPointer<vtkInEdgeIterator>::New();
	reebGraph->GetInEdges(nodeId, inEdges);
	while(inEdges->HasNext()){
		adjacentVertices.push_back(inEdges->Next().Source);
	}
	vtkSmartPointer<vtkOutEdgeIterator> outEdges = vtkSmartPointer<vtkOutEdgeIterator>::New();
	reebGraph->GetOutEdges(nodeId, outEdges);
	while(outEdges->HasNext()){
		adjacentVertices.push_back(outEdges->Next().Target);
	}
	int maximumNode = -1;
	for(int i = 0; i < adjacentVertices.size(); i++){
		int adjacentnodeId = adjacentVertices[i];
		double adjacentNodeRange = field->GetComponent(adjacentnodeId, 0);
		if(find(path.begin(), path.end(), adjacentnodeId) == path.end() && adjacentNodeRange > minimumRangeValue){
			int currentMaximumNode = GetMaximumNode(adjacentnodeId, minimumRangeValue, path);
			if(maximumNode == -1){
				maximumNode = currentMaximumNode;
			}
			else if(currentMaximumNode != -1){
				if(field->GetComponent(currentMaximumNode, 0) > field->GetComponent(maximumNode, 0) ){
					maximumNode = currentMaximumNode;
				}		
			}
		}
	}
	return maximumNode;
}


//Starting from "nodeId", traverse all the nodes of the Reeb graph having their range values between "min_range" and "max_range"
void vtkPersistenceDiagramOfReebGraph::VisitNodesWithinASpecifiedRange(int nodeId, double min_range, double max_range, vector<bool> &visited){
	visited[nodeId] = true;
	vector<int> adjacentVertices;
	vtkSmartPointer<vtkDoubleArray> field;
	field = vtkDoubleArray::SafeDownCast(reebGraph->GetVertexData()->GetArray(fieldName));
	vtkSmartPointer<vtkInEdgeIterator> inEdges = vtkSmartPointer<vtkInEdgeIterator>::New();
	reebGraph->GetInEdges(nodeId, inEdges);
	while(inEdges->HasNext()){
		adjacentVertices.push_back(inEdges->Next().Source);
	}
	vtkSmartPointer<vtkOutEdgeIterator> outEdges = vtkSmartPointer<vtkOutEdgeIterator>::New();
	reebGraph->GetOutEdges(nodeId, outEdges);
	while(outEdges->HasNext()){
		adjacentVertices.push_back(outEdges->Next().Target);
	}
	for(int i = 0; i < adjacentVertices.size(); i++){
		int adjacentNodeId = adjacentVertices[i];
		double adjacentNodeRange = field->GetComponent(adjacentNodeId, 0);
		if(adjacentNodeRange > min_range && adjacentNodeRange < max_range){
			if(visited[adjacentNodeId] == false){
				VisitNodesWithinASpecifiedRange(adjacentNodeId, min_range, max_range, visited);
			}
		}
	}
}

//Constructs the four types of persistence diagrams
void vtkPersistenceDiagramOfReebGraph::ConstructPersistenceDiagram(){
	vtkSmartPointer<vtkDoubleArray> field;
	field = vtkDoubleArray::SafeDownCast(reebGraph->GetVertexData()->GetArray(fieldName));
	int nrVertices = reebGraph->GetNumberOfVertices();
	int globalMinimum = -1, globalMaximum = -1;
	for(int i = 0; i < nrVertices; i++){
		int inDegree = reebGraph->GetInDegree(i);
		int outDegree = reebGraph->GetOutDegree(i);
		if(outDegree == 2){//Up-fork
			vtkSmartPointer<vtkOutEdgeIterator> edges = vtkSmartPointer<vtkOutEdgeIterator>::New();
			reebGraph->GetOutEdges(i, edges);
			int adjacentnodeId1 = edges->Next().Target;
			int adjacentnodeId2 = edges->Next().Target;
			vector<int> path1;
			int maximumNodeOfComponent1 = GetMaximumNode(adjacentnodeId1, field->GetComponent(i, 0), path1);
			vector<int> path2;
			int maximumNodeOfComponent2 = GetMaximumNode(adjacentnodeId2, field->GetComponent(i, 0), path2);
			if(maximumNodeOfComponent1 == maximumNodeOfComponent2){//Essential up-fork
				//repetition - pairing is not required (taken care while processing the essential down-fork)
			}
			else{//Ordinary up-fork - pair it with the appropriate maxmimum node
				if(field->GetComponent(maximumNodeOfComponent1, 0) > field->GetComponent(maximumNodeOfComponent2, 0)){
					persistenceDiagramSuperLevelSetFiltration.push_back(make_pair(field->GetComponent(i, 0), field->GetComponent(maximumNodeOfComponent2, 0)));
				}
				else{
					persistenceDiagramSuperLevelSetFiltration.push_back(make_pair(field->GetComponent(i, 0), field->GetComponent(maximumNodeOfComponent1, 0)));
				}
			}
		}
		else if(inDegree == 2 ){
			vtkSmartPointer<vtkInEdgeIterator> edges = vtkSmartPointer<vtkInEdgeIterator>::New();
			reebGraph->GetInEdges(i, edges);
			int adjacentnodeId1 = edges->Next().Source;
			int adjacentnodeId2 = edges->Next().Source;
			vector<int> path1;			
			int minimumNodeOfComponent1 = GetMinimumNode(adjacentnodeId1, field->GetComponent(i, 0), path1);
			vector<int> path2;
			int minimumNodeOfComponent2 = GetMinimumNode(adjacentnodeId2, field->GetComponent(i, 0), path2);
			if(minimumNodeOfComponent1 == minimumNodeOfComponent2){//Essential down-fork
				double cur_max_up_fork_value = INT_MIN;//for storing the up-fork corresponding to the smallest cycle
				for(int v_id = 0; v_id < reebGraph->GetNumberOfVertices(); v_id++){
					if(reebGraph->GetOutDegree(v_id) == 2){
						vector<bool> visited_array;
						for(int t_i = 0; t_i < reebGraph->GetNumberOfVertices(); t_i++){
							visited_array.push_back(false);
						}
						VisitNodesWithinASpecifiedRange(v_id, field->GetComponent(v_id,0), field->GetComponent(i, 0),visited_array);
						if(visited_array[adjacentnodeId1] == true && visited_array[adjacentnodeId2] == true && cur_max_up_fork_value < field->GetComponent(v_id,0)){
							cur_max_up_fork_value = field->GetComponent(v_id,0);
						}
					}
				}
				extendedPersistenceDiagram_dimension_1.push_back(make_pair(field->GetComponent(i, 0),cur_max_up_fork_value));
			}
			else{//Ordinary down-fork - pair it with the appropriate minmimum node
				if(field->GetComponent(minimumNodeOfComponent1, 0) < field->GetComponent(minimumNodeOfComponent2, 0)){
					persistenceDiagramSublevelSetFiltration.push_back(make_pair(field->GetComponent(minimumNodeOfComponent2, 0), field->GetComponent(i, 0)));
				}

				else{
					persistenceDiagramSublevelSetFiltration.push_back(make_pair(field->GetComponent(minimumNodeOfComponent1, 0), field->GetComponent(i, 0)));
				}
			}
		}
		else if(inDegree == 0){//finding the global minimum
			if(globalMinimum == -1){
				globalMinimum = i;
			}
			else if(field->GetComponent(i, 0) < field->GetComponent(globalMinimum, 0)){
				globalMinimum = i;
			}
		}
		else if(outDegree == 0){//finding the global maximum
			if(globalMaximum == -1){
				globalMaximum = i;
			}
			else if(field->GetComponent(i, 0) > field->GetComponent(globalMaximum, 0)){
				globalMaximum = i;
			}
		}
	}
	if(nrVertices > 1){
		//Add the point "(globalMinimum, globalMaximum)"
		persistenceDiagramSublevelSetFiltration.push_back(make_pair(field->GetComponent(globalMinimum, 0),field->GetComponent(globalMaximum, 0)));
		persistenceDiagramSuperLevelSetFiltration.push_back(make_pair(field->GetComponent(globalMinimum, 0), field->GetComponent(globalMaximum, 0)));
		extendedPersistenceDiagram_dimension_0.push_back(make_pair(field->GetComponent(globalMinimum, 0), field->GetComponent(globalMaximum, 0)));
	}
	else{
		//Add the point "(globalMinimum, globalMaximum)" (in this case globalMinimum = globalMaximum)
		persistenceDiagramSublevelSetFiltration.push_back(make_pair(field->GetComponent(globalMinimum, 0),field->GetComponent(globalMinimum, 0)));
		persistenceDiagramSuperLevelSetFiltration.push_back(make_pair(field->GetComponent(globalMinimum, 0), field->GetComponent(globalMinimum, 0)));
		extendedPersistenceDiagram_dimension_0.push_back(make_pair(field->GetComponent(globalMinimum, 0), field->GetComponent(globalMinimum, 0)));//Global minimum and maximum is the same.
	}
}
