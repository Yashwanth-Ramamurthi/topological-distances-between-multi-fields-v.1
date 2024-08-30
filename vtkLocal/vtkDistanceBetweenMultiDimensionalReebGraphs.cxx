/*=========================================================================

 *	File: vtkDistanceBetweenMultiDimensionalReebGraphs.cxx
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
//Reference:
//[1] Yashwanth Ramamurthi and Amit Chattopadhyay. 2023. Topological Shape Matching using Multi-Dimensional Reeb Graphs. In Proceedings of the Thirteenth Indian Conference on Computer Vision, Graphics and Image Processing (ICVGIP '22). Association for Computing Machinery, New York, NY, USA, Article 5, 1–10.

#include "vtkObjectFactory.h"
#include "vtkDoubleArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkInEdgeIterator.h"
#include "vtkOutEdgeIterator.h"
#include "vtkAdjacentVertexIterator.h"
#include "vtkIdTypeArray.h"
#include "vtkDistanceBetweenMultiDimensionalReebGraphs.h"
vtkStandardNewMacro(vtkDistanceBetweenMultiDimensionalReebGraphs);

#define VTK_CREATE(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

vtkDistanceBetweenMultiDimensionalReebGraphs::vtkDistanceBetweenMultiDimensionalReebGraphs(){
    this->Fields = std::vector<char*>();
}

vtkDistanceBetweenMultiDimensionalReebGraphs::~vtkDistanceBetweenMultiDimensionalReebGraphs(){
}

// Add the names of the component scalar fields of the bivariate fields (eg., f1,f2).
void vtkDistanceBetweenMultiDimensionalReebGraphs::AddField(char *fieldName){
	char *f = new char[strlen(fieldName) + 1];
	strcpy(f, fieldName);
	this->Fields.push_back(f);
}

//Segregation of the persistence diagrams of the component Reeb graphs of the MDRG of a multi-field f = (f1,f2,...,fn), into sets of the form S_{f_i}^{c}
//persistenceDiagramMap[i][c] stores the list of persistence diagrams S_{f_i}^{c} (see Eqution (12) in [1])
void vtkDistanceBetweenMultiDimensionalReebGraphs::ConstructPersistenceDiagramMap(vtkMultiDimensionalReebGraph *mdrg, int graphId, int fieldId, vector< map<double, vector<vector<pair<double,double> > > > > &persistenceDiagramMap, int type){
	
	vtkReebGraph* reebGraph = mdrg->GetTreeNode(graphId)->GetReebGraph();
	char *fieldName = mdrg->GetTreeNode(graphId)->GetFieldName();
	vtkSmartPointer<vtkDoubleArray> fieldArray;
	fieldArray = vtkDoubleArray::SafeDownCast(reebGraph->GetVertexData()->GetArray(fieldName));
	
	vtkSmartPointer<vtkAdjacentVertexIterator> adj = vtkAdjacentVertexIterator::New();
	mdrg->GetMultiDimensionalReebGraph()->GetAdjacentVertices(graphId, adj);
	
	int i = 0;
	while(adj->HasNext()){
		int childGraphId = adj->Next();
		ConstructPersistenceDiagramMap(mdrg, childGraphId, fieldId+1, persistenceDiagramMap, type);
		vtkReebGraph* childReebGraph = mdrg->GetTreeNode(childGraphId)->GetReebGraph();
		char *childReebGraphFieldName = mdrg->GetTreeNode(childGraphId)->GetFieldName();
		//Handle degeneracy in the Reeb graph
		vtkHandleDegeneracyInReebGraph *handleDegeneracy = vtkHandleDegeneracyInReebGraph::New();
		handleDegeneracy->SetParameters(childReebGraph, childReebGraphFieldName);
		handleDegeneracy->HandleDegeneracy();
	
		//Construct the persistence diagram
		vtkPersistenceDiagramOfReebGraph* persistenceDiagramOfChildReebGraph = vtkPersistenceDiagramOfReebGraph::New();
		persistenceDiagramOfChildReebGraph->fieldName = childReebGraphFieldName;
		persistenceDiagramOfChildReebGraph->reebGraph = childReebGraph;
		persistenceDiagramOfChildReebGraph->ConstructPersistenceDiagram();
		vector<pair<double, double> > persistenceDiagram = persistenceDiagramOfChildReebGraph->GetPersistenceDiagram(type);
		if(persistenceDiagramMap[fieldId + 1].find(fieldArray->GetValue(i)) == persistenceDiagramMap[fieldId + 1].end()){
			vector<vector<pair<double, double> > > persistenceDiagrams;
			persistenceDiagrams.push_back(persistenceDiagram);
			persistenceDiagramMap[fieldId + 1][fieldArray->GetValue(i)] = persistenceDiagrams;
		}
		else{
			persistenceDiagramMap[fieldId + 1][fieldArray->GetValue(i)].push_back(persistenceDiagram);
		}

		i++;
	}
	
	if(fieldId == 0){
		//Handle degeneracy in the Reeb graph
		vtkHandleDegeneracyInReebGraph *handleDegeneracy = vtkHandleDegeneracyInReebGraph::New();
		handleDegeneracy->SetParameters(reebGraph, fieldName);
		handleDegeneracy->HandleDegeneracy();

		//Construct the persistence diagram
		vtkPersistenceDiagramOfReebGraph* persistenceDiagramOfReebGraph = vtkPersistenceDiagramOfReebGraph::New();
		persistenceDiagramOfReebGraph->fieldName = fieldName;
		persistenceDiagramOfReebGraph->reebGraph = reebGraph;
		persistenceDiagramOfReebGraph->ConstructPersistenceDiagram();
		vector<pair<double, double> > persistenceDiagram = persistenceDiagramOfReebGraph->GetPersistenceDiagram(type);
		vector<vector<pair<double, double> > > persistenceDiagrams;
		persistenceDiagrams.push_back(persistenceDiagram);
		persistenceDiagramMap[0][0] = persistenceDiagrams;	
	}
}

//Computing the minimum-cost matching for the computation of the bottleneck distance (based on the Hungarian algorithm)
double vtkDistanceBetweenMultiDimensionalReebGraphs::MinCostMatching(const vector<vector<double> > &cost, vector<int> &Lmate, vector<int> &Rmate) {
  int n = int(cost.size());

  // construct dual feasible solution
  vector<double> u(n);
  vector<double> v(n);
  for (int i = 0; i < n; i++) {
    u[i] = cost[i][0];
    for (int j = 1; j < n; j++) u[i] = min(u[i], cost[i][j]);
  }
  for (int j = 0; j < n; j++) {
    v[j] = cost[0][j] - u[0];
    for (int i = 1; i < n; i++) v[j] = min(v[j], cost[i][j] - u[i]);
  }
  
  // construct primal solution satisfying complementary slackness
  Lmate = vector<int>(n, -1);
  Rmate = vector<int>(n, -1);
  int mated = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (Rmate[j] != -1) continue;
      if (fabs(cost[i][j] - u[i] - v[j]) < 1e-10) {
        Lmate[i] = j;
        Rmate[j] = i;
        mated++;
        break;
      }
    }
  }
  
  vector<double> dist(n);
  vector<int> dad(n);
  vector<int> seen(n);
  
  // repeat until primal solution is feasible
  while (mated < n) {
    
    // find an unmatched left node
    int s = 0;
    while (Lmate[s] != -1) s++;
    
    // initialize Dijkstra
    fill(dad.begin(), dad.end(), -1);
    fill(seen.begin(), seen.end(), 0);
    for (int k = 0; k < n; k++) 
      dist[k] = cost[s][k] - u[s] - v[k];
    
    int j = 0;
    while (true) {
      
      // find closest
      j = -1;
      for (int k = 0; k < n; k++) {
        if (seen[k]) continue;
        if (j == -1 || dist[k] < dist[j]) j = k;
      }
      seen[j] = 1;
      
      // termination condition
      if (Rmate[j] == -1) break;
      
      // relax neighbors
      const int i = Rmate[j];
      for (int k = 0; k < n; k++) {
        if (seen[k]) continue;
        const double new_dist = dist[j] + cost[i][k] - u[i] - v[k];
        if (dist[k] > new_dist) {
          dist[k] = new_dist;
          dad[k] = j;
        }
      }
    }
    
    // update dual variables
    for (int k = 0; k < n; k++) {
      if (k == j || !seen[k]) continue;
      const int i = Rmate[k];
      v[k] += dist[k] - dist[j];
      u[i] -= dist[k] - dist[j];
    }
    u[s] += dist[j];
    
    while (dad[j] >= 0) {
      const int d = dad[j];
      Rmate[j] = Rmate[d];
      Lmate[Rmate[j]] = j;
      j = d;
    }
    Rmate[j] = s;
    Lmate[s] = j;
    
    mated++;
  }
  
  double value = 0;
  for (int i = 0; i < n; i++)
    value += cost[i][Lmate[i]];
  
  return value;
}

//Compute the bottleneck Distance between two persistence diagrams
double vtkDistanceBetweenMultiDimensionalReebGraphs::BottleneckDistance(vector<pair<double, double> > persistenceDiagramOfReebGraph1, vector<pair<double, double> > persistenceDiagramOfReebGraph2){

	int numberOfPoints = persistenceDiagramOfReebGraph1.size() + persistenceDiagramOfReebGraph2.size();
	if(numberOfPoints == 0){
		return 0;
	}
	std::vector<double> v(numberOfPoints, INT_MAX);//points + diagonal points
	std::vector<std::vector<double> > costMatrix(numberOfPoints, v);
	map<pair<double, double>, vector<vector<int> > >::iterator itr1;
	int i = 0;

  	for (int i = 0; i < persistenceDiagramOfReebGraph1.size(); i++) {
        	pair<double, double> point1 = persistenceDiagramOfReebGraph1[i];
        	int j = 0;
	  	for (int j = 0; j < persistenceDiagramOfReebGraph2.size(); j++) {
	        	pair<double, double> point2 = persistenceDiagramOfReebGraph2[j];
			costMatrix[i][j] = max(abs(point2.first - point1.first), abs(point2.second - point1.second));
        	}
	}
	i = 0;
  	for (int i = 0; i < persistenceDiagramOfReebGraph1.size(); i++) {
        	pair<double, double> point = persistenceDiagramOfReebGraph1[i];
		costMatrix[i][persistenceDiagramOfReebGraph2.size() + i] = ((abs(point.second - point.first))/2.0);
	}
	i = 0;
        map<pair<double, double>, vector<vector<int> > >::iterator itr2;
  	for (int i = 0; i < persistenceDiagramOfReebGraph2.size(); i++) {
        	pair<double, double> point = persistenceDiagramOfReebGraph2[i];
		costMatrix[persistenceDiagramOfReebGraph1.size() + i][i] = ((abs(point.second - point.first))/2.0);
	}

	for(int i = persistenceDiagramOfReebGraph1.size(); i < numberOfPoints; i++){
		for(int j = persistenceDiagramOfReebGraph2.size(); j < numberOfPoints; j++){
			costMatrix[i][j]= 0;//Matching diagonal points should take zero cost
		}
	}
	vector<double> costs;
	for(int i = 0; i < costMatrix.size();i++){
		for(int j = 0; j < costMatrix[i].size(); j++){
			costs.push_back(costMatrix[i][j]);
		}
	}
	sort(costs.begin(), costs.end());
	int l_index = 0;
	int r_index = costs.size()-1;
	while(r_index > l_index){
		int m_index = floor((l_index + r_index)/2);
		double thresh = costs[m_index];
		vector <vector<double> > modifiedCostMatrix = costMatrix;
		for(int i = 0; i < modifiedCostMatrix.size();i++){
			for(int j = 0; j < modifiedCostMatrix[i].size(); j++){
				if(modifiedCostMatrix[i][j] <= thresh){
					modifiedCostMatrix[i][j] = 0;
				}
				else{
					modifiedCostMatrix[i][j] = 1;
				}
			}
		}
		vector<int> Lmate;
		vector<int> Rmate;
		double cost = MinCostMatching(modifiedCostMatrix, Lmate, Rmate);
		
		if(cost == 0){
			r_index = m_index;
		}
		else{
			if(l_index == m_index){
				l_index++;
			}
			else{
				l_index = m_index;
			}
		}
	}
	double bottleneckDistance = costs[l_index];
	return bottleneckDistance;	
}

//Compute the bottleneck distance between the persistence diagram "persistenceDiagramOfReebGraph1" and a persistence diagram consisting of infinitely many diagonal points (Each point in "persistence diagram" is matched to the closest diagonal point)
double vtkDistanceBetweenMultiDimensionalReebGraphs::BottleneckDistance(vector<pair<double, double> > persistenceDiagram){
	double distance = 0;
	for(int i = 0; i < persistenceDiagram.size(); i++){
		double val = abs(persistenceDiagram[i].second - persistenceDiagram[i].first)/2.0;
		distance = max(distance, val);
	}
	return distance;
}

//Compute the sum of the distances between the sets of persistence diagrams in "persistenceDiagrams1" and "persistenceDiagrams2"
//persistenceDiagrams1: consists of the persistence diagrams corresponding to the Reeb graphs in S_{f_{i+1}}^{c} (see Eqution (12) in [1])
//persistenceDiagrams2: consists of the persistence diagrams corresponding to the Reeb graphs in S_{g_{i+1}}^{c} (see Eqution (12) in [1])
double vtkDistanceBetweenMultiDimensionalReebGraphs:: DistanceBetweenSetsOfPersistenceDiagrams(vector<vector<pair<double,double> > > &persistenceDiagrams1, vector<vector<pair<double,double> > > &persistenceDiagrams2){
	int n1 = persistenceDiagrams1.size();
	int n2 = persistenceDiagrams2.size();
	int size = n1+n2;
	std::vector<double> v(size, INT_MAX);
	std::vector<std::vector<double> > costMatrix(size, v);
	for(int i = 0; i < n1; i++){
		for(int j = 0; j < n2; j++){
			costMatrix[i][j] = BottleneckDistance(persistenceDiagrams1[i], persistenceDiagrams2[j]);
		}
	}
	for(int i = 0 ; i < n1; i++){
		double val = BottleneckDistance(persistenceDiagrams1[i]);
		costMatrix[i][n2+i] =val;
	}
	for(int j = 0 ; j < n2; j++){
		double val = BottleneckDistance(persistenceDiagrams2[j]);
		costMatrix[n1+j][j] =val;
	}
	for(int i = n1 ; i < size; i++){
		for(int j = n2; j < size; j++){
			costMatrix[i][j]= 0;//Matching diagonal points should take zero cost
		}
	}
	vector<double> costs;
	for(int i = 0; i < costMatrix.size();i++){
		for(int j = 0; j < costMatrix[i].size(); j++){
			costs.push_back(costMatrix[i][j]);
		}
	}
	sort(costs.begin(), costs.end());
	int l_index = 0;
	int r_index = costs.size()-1;
	while(r_index > l_index){
		int m_index = floor((l_index + r_index)/2);
		double thresh = costs[m_index];
		vector <vector<double> > modifiedCostMatrix = costMatrix;
		for(int i = 0; i < modifiedCostMatrix.size();i++){
			for(int j = 0; j < modifiedCostMatrix[i].size(); j++){
				if(modifiedCostMatrix[i][j] <= thresh){
					modifiedCostMatrix[i][j] = 0;
				}
				else{
					modifiedCostMatrix[i][j] = 1;
				}
			}
		}
		vector<int> Lmate;
		vector<int> Rmate;
		double cost = MinCostMatching(modifiedCostMatrix, Lmate, Rmate);
		
		if(cost == 0){
			r_index = m_index;
		}
		else{
			if(l_index == m_index){
				l_index++;
			}
			else{
				l_index = m_index;
			}
		}
	}
	return costs[l_index];	
}

//Compute the distance between the MDRGs corresponding to the JCNs "JCN1" and "JCN2" based on the persistence diagrams of a particular "type"
//type 0 => 0-dimensional Persistence Diagram corresponding to the sublevelset filtration 
//type 1 => 0-dimensional Persistence Diagram corresponding to the superlevelset filtration
//type 2 => 0-dimensional Extended Persistence Diagram
//type 3 => 1-dimensional Extended Persistence Diagram
//Note: Distances computed for types 0, 1, and 3 need to be multiplied by the weights w_0, w_1, and w_2, respectively, and the resulting values need to be be summed to obtain the distance d_T (see Equations 10 and 11 in [1])
double vtkDistanceBetweenMultiDimensionalReebGraphs::ComputeDistance(vtkGraph * JCN1, vtkGraph * JCN2, int type){
	persistenceDiagramMap1 = *(new vector< map<double, vector<vector<pair<double,double> > > > >(this->Fields.size()));
	persistenceDiagramMap2 = *(new vector< map<double, vector<vector<pair<double,double> > > > >(this->Fields.size()));
	int nrVerticesJCN1 = JCN1->GetNumberOfVertices();
	VTK_CREATE(vtkIdTypeArray, nodeIdsJCN1);
	nodeIdsJCN1->SetName("jcnids");
	nodeIdsJCN1->SetNumberOfValues(nrVerticesJCN1);
	for (vtkIdType i = 0; i < nrVerticesJCN1; i++)
	{
	     nodeIdsJCN1->SetValue(i, i);//copying the origical JCN-ids  
    	}
 	JCN1->GetVertexData()->AddArray(nodeIdsJCN1);//Adding JCN-ids in Input JCN 1
	int nrVerticesJCN2 = JCN2->GetNumberOfVertices();
	VTK_CREATE(vtkIdTypeArray, nodeIdsJCN2);
	nodeIdsJCN2->SetName("jcnids");
	nodeIdsJCN2->SetNumberOfValues(nrVerticesJCN2);
	for (vtkIdType i = 0; i < nrVerticesJCN2; i++)
	{
	     nodeIdsJCN2->SetValue(i, i);//copying the origical JCN-ids  
    	}
 	JCN2->GetVertexData()->AddArray(nodeIdsJCN2);//Adding JCN-ids in Input JCN 2

	vtkMultiDimensionalReebGraph* mdrg1 = vtkMultiDimensionalReebGraph::New();
	vtkMultiDimensionalReebGraph* mdrg2 = vtkMultiDimensionalReebGraph::New();	
	for(int f = 0; f < this->Fields.size(); f++){
		mdrg1->AddField(this->Fields[f]);
		mdrg2->AddField(this->Fields[f]);
	}
	mdrg1->CreateMultiDimensionalReebGraph(JCN1);
	mdrg2->CreateMultiDimensionalReebGraph(JCN2);
	ConstructPersistenceDiagramMap(mdrg1, 0, 0, persistenceDiagramMap1, type);
	ConstructPersistenceDiagramMap(mdrg2, 0, 0, persistenceDiagramMap2, type);
	double distance_all_fields = 0;
	for(int f = 0; f < this->Fields.size(); f++){
		std::set<double> fieldValues;
		for(map<double, vector<vector<pair<double,double> > > >::iterator it = persistenceDiagramMap1[f].begin(); it != persistenceDiagramMap1[f].end(); ++it){
		 	 fieldValues.insert(it->first);
		}
		for(map<double, vector<vector<pair<double,double> > > >::iterator it = persistenceDiagramMap2[f].begin(); it != persistenceDiagramMap2[f].end(); ++it){
		 	 fieldValues.insert(it->first);
		}
		double distance_currentField = 0;
		for(set<double>::iterator it = fieldValues.begin(); it != fieldValues.end(); ++it) {
			vector<vector<pair<double,double> > > persistenceDiagramList1(1), persistenceDiagramList2(1);
			if(persistenceDiagramMap1[f].find(*it) != persistenceDiagramMap1[f].end()){
				persistenceDiagramList1 = persistenceDiagramMap1[f][*it];
			}
			if(persistenceDiagramMap2[f].find(*it) != persistenceDiagramMap2[f].end()){
				persistenceDiagramList2 = persistenceDiagramMap2[f][*it];
			}
			double distance;
			distance = DistanceBetweenSetsOfPersistenceDiagrams(persistenceDiagramList1, persistenceDiagramList2);
			distance_currentField = distance_currentField + distance;
		}
		distance_all_fields = distance_all_fields + (distance_currentField/fieldValues.size());
	}
	return distance_all_fields;
}
