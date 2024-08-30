/*=========================================================================

 *	File: vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs.cxx
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
//Reference:
//[1] Yashwanth Ramamurthi and Amit Chattopadhyay, "A Topological Distance Between Multi-Fields Based on Multi-Dimensional Persistence Diagrams," in IEEE Transactions on Visualization and Computer
//Graphics, vol. 30, no. 9, pp. 5939-5952, Sept. 2024

#include "vtkObjectFactory.h"
#include "vtkDoubleArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkInEdgeIterator.h"
#include "vtkOutEdgeIterator.h"
#include "vtkAdjacentVertexIterator.h"
#include "vtkIdTypeArray.h"
#include "vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs.h"
vtkStandardNewMacro(vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs);

#define VTK_CREATE(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs(){
    this->Fields = std::vector<char*>();
}

vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::~vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs(){
}

// Add the names of the component scalar fields of the bivariate fields (eg., f1,f2).
void vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::AddField(char *fieldName){
	char *f = new char[strlen(fieldName) + 1];
	strcpy(f, fieldName);
	this->Fields.push_back(f);
}

// Set the parameter q for computing the Wasserstein distance (see Equations 10 and 11 in [1])
void vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::SetParameterQ(double q){
	this->q = q;
}

//Convert a double (floating point) value to a string
string vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::DoubleToString(double number){
    std::ostringstream out;
    out<<std::setprecision(20)<<number;
    return out.str();
}


//Computes the Multi-Dimensional Persistence Diagram PD_{f} corresponding to an MDRG MR_{f}. 
// The points in the MDPD are stored in groups of the form PD_p^{x_1}(f)//
// For a real-number v, multiDimensionalPersistenceDiagramMap[v] stores a list of {PD_p^{x_1}(f): \overline{f_1}(p) = val} (see Equation 4 in [1] for the definition of PD_p^{x_1}(f))
// The segregation of the points in the MDPD into groups helps in the checking of criteria (C1)-(C3) while computing the Wasserstein distance between two MDPDs
void vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::ConstructPersistenceDiagramMap(vtkMultiDimensionalReebGraph *mdrg, int graphId, int fieldId, map<string, vector<vector<vector<double> > > > &multiDimensionalPersistenceDiagramMap, vector<int> &types, vector<vector<double>> combinedParentPersistenceDiagram, string fieldValue, double parentFieldValue){
	vtkReebGraph* reebGraph = mdrg->GetTreeNode(graphId)->GetReebGraph();
	char *fieldName = mdrg->GetTreeNode(graphId)->GetFieldName();
	vtkSmartPointer<vtkDoubleArray> fieldArray;
	fieldArray = vtkDoubleArray::SafeDownCast(reebGraph->GetVertexData()->GetArray(fieldName));
	vtkSmartPointer<vtkIdTypeArray> jcnNodeIdArray;
	jcnNodeIdArray = vtkIdTypeArray::SafeDownCast(reebGraph->GetVertexData()->GetArray("jcnids"));
	int numberOfVerticesBeforeRemovingDegenaracy = reebGraph->GetNumberOfVertices();
	vtkHandleDegeneracyInReebGraph *handleDegeneracy = vtkHandleDegeneracyInReebGraph::New();
	handleDegeneracy->SetParameters(reebGraph, fieldName);
	handleDegeneracy->HandleDegeneracy();
	vtkPersistenceDiagramOfReebGraph* persistenceDiagramOfReebGraph = vtkPersistenceDiagramOfReebGraph::New();
	persistenceDiagramOfReebGraph->fieldName = fieldName;
	persistenceDiagramOfReebGraph->reebGraph = reebGraph;
	persistenceDiagramOfReebGraph->ConstructPersistenceDiagram();
	vector<pair<double, double> > persistenceDiagram;
	for(int t = 0; t < types.size(); t++){
		vector<pair<double, double> > persistenceDiagram_t = persistenceDiagramOfReebGraph->GetPersistenceDiagram(types[t]);
		persistenceDiagram.insert(persistenceDiagram.end(), persistenceDiagram_t.begin(), persistenceDiagram_t.end());
	}
	persistenceDiagram.push_back(make_pair(0.0,0.0));// Adding a diagonal point to make sure that if this diagram has zero points it does not affect the computation of the Wasserstein distance in the next dimension.
	vector<vector<double> > combinedPersistenceDiagram;
	if(fieldId == 0){
		for(int k = 0; k < persistenceDiagram.size(); k++){
			vector<double> point;
			point.push_back(persistenceDiagram[k].first);
			point.push_back(persistenceDiagram[k].second);
			combinedPersistenceDiagram.push_back(point);
		}
	}
	else{
		for(int j = 0; j < combinedParentPersistenceDiagram.size(); j++){
			for(int k = 0; k < persistenceDiagram.size(); k++){
				vector<double> point_combinedPersistenceDiagram = combinedParentPersistenceDiagram[j];
				double parent_birth = point_combinedPersistenceDiagram[point_combinedPersistenceDiagram.size()-2];
				double parent_death = point_combinedPersistenceDiagram[point_combinedPersistenceDiagram.size()-1];
				point_combinedPersistenceDiagram.push_back(persistenceDiagram[k].first);
				point_combinedPersistenceDiagram.push_back(persistenceDiagram[k].second);
				if(min(parent_birth,parent_death) <= parentFieldValue && parentFieldValue <= max(parent_birth,parent_death)){
					combinedPersistenceDiagram.push_back(point_combinedPersistenceDiagram);
				}
			}
		}
	}
	if(fieldId == (this->Fields.size() - 1)){
		if(multiDimensionalPersistenceDiagramMap.find(fieldValue) == multiDimensionalPersistenceDiagramMap.end()){
			vector<vector<vector<double > > > combinedPersistenceDiagrams;
			combinedPersistenceDiagrams.push_back(combinedPersistenceDiagram);
			multiDimensionalPersistenceDiagramMap[fieldValue] = combinedPersistenceDiagrams;
		}
		else{
			multiDimensionalPersistenceDiagramMap[fieldValue].push_back(combinedPersistenceDiagram);
		}
	}
	
	vtkSmartPointer<vtkAdjacentVertexIterator> adj = vtkAdjacentVertexIterator::New();
	mdrg->GetMultiDimensionalReebGraph()->GetAdjacentVertices(graphId, adj);
	fieldArray = vtkDoubleArray::SafeDownCast(reebGraph->GetVertexData()->GetArray(fieldName));
	if(adj->HasNext()){
		for(int i = 0; i <  numberOfVerticesBeforeRemovingDegenaracy; i++){
			int childGraphId = adj->Next();
			string s = fieldValue + " " +  DoubleToString(fieldArray->GetValue(i));
			ConstructPersistenceDiagramMap(mdrg, childGraphId, fieldId+1, multiDimensionalPersistenceDiagramMap, types, combinedPersistenceDiagram, s, fieldArray->GetValue(i));
		}
	}
}

//Computing the minimum-cost matching for the computation of the Wasserstein distance (based on the Hungarian algorithm)
double vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::MinCostMatching(const vector<vector<double> > &cost, vector<int> &Lmate, vector<int> &Rmate) {
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




//Compute the Wasserstein distance between sets of points in the MDPD grouped into different categories based on the matching criteria (C1)-(C3)
double vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::WassersteinDistanceBetweenSetsOfPeristenceDiagrams(vector<vector<vector<double> > > &persistenceDiagrams1, vector<vector<vector<double> > > &persistenceDiagrams2){
	
	int numberOfDiagrams = persistenceDiagrams1.size() + persistenceDiagrams2.size();
	if(numberOfDiagrams == 0){
		return 0;
	}
	std::vector<double> v(numberOfDiagrams, INT_MAX);//points + diagonal points (persistence diagrams)
	std::vector<std::vector<double> > costMatrix(numberOfDiagrams, v);
  	for (int i = 0; i < persistenceDiagrams1.size(); i++) {
        	vector<vector<double> > &persistenceDiagram1 = persistenceDiagrams1[i];
	  	for (int j = 0; j < persistenceDiagrams2.size(); j++) {
        		vector<vector<double> > &persistenceDiagram2 = persistenceDiagrams2[j];
			costMatrix[i][j] = WassersteinDistance(persistenceDiagram1,persistenceDiagram2);
        	}
	}
  	for (int i = 0; i < persistenceDiagrams1.size(); i++) {
  	        vector<vector<double> > &persistenceDiagram = persistenceDiagrams1[i];
  	        costMatrix[i][persistenceDiagrams2.size() + i] = 0;
  	        for(int c = 0; c < persistenceDiagram.size(); c++){
			double cost = 0;
			for(int p = 0; p < persistenceDiagram[c].size(); p = p + 2){
				cost = max(cost, (abs(persistenceDiagram[c][p+1]- persistenceDiagram[c][p])/2.0));
			}
			costMatrix[i][persistenceDiagrams2.size() + i] += pow(cost, q);
  	        }
	}
  	
  	for (int i = 0; i < persistenceDiagrams2.size(); i++) {
  	        vector<vector<double> > &persistenceDiagram = persistenceDiagrams2[i];
  	        costMatrix[persistenceDiagrams1.size() + i][i] = 0;
  	        for(int c = 0; c < persistenceDiagram.size(); c++){
			double cost = 0;
			for(int p = 0; p < persistenceDiagram[c].size(); p = p + 2){
				cost = max(cost, (abs(persistenceDiagram[c][p+1]- persistenceDiagram[c][p])/2.0));
	  	        }
			costMatrix[persistenceDiagrams1.size() + i][i] += pow(cost, q);
		}
	}

	for(int i = persistenceDiagrams1.size(); i < numberOfDiagrams; i++){
		for(int j = persistenceDiagrams2.size(); j < numberOfDiagrams; j++){
			costMatrix[i][j]= 0;//Matching diagonal points (persistence diagrams) should take zero cost
		}
	}
	vector<int> Lmate;
	vector<int> Rmate;
	double cost = MinCostMatching(costMatrix, Lmate, Rmate);
	return cost;	
}

//Compute the Wasserstein Distance between two persistence diagrams
double vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::WassersteinDistance(vector<vector<double> > &persistenceDiagram1, vector<vector<double> > &persistenceDiagram2){
	int numberOfPoints = persistenceDiagram1.size() + persistenceDiagram2.size();
	if(numberOfPoints == 0){
		return 0;
	}
	std::vector<double> v(numberOfPoints, INT_MAX);//points + diagonal points
	std::vector<std::vector<double> > costMatrix(numberOfPoints, v);
	map<pair<double, double>, vector<vector<int> > >::iterator itr1;
	int i = 0;
  	for (int i = 0; i < persistenceDiagram1.size(); i++) {
        	vector<double> &point1 = persistenceDiagram1[i];
        	int j = 0;
	  	for (int j = 0; j < persistenceDiagram2.size(); j++) {
	        	vector<double> &point2 = persistenceDiagram2[j];
			costMatrix[i][j] = 0;
			for(int k = 0; k < point1.size(); k++){
				costMatrix[i][j] = max(costMatrix[i][j], abs(point2[k] - point1[k]));
			}
			costMatrix[i][j] = pow(costMatrix[i][j], q);
        	}
	}
	i = 0;
  	for (int i = 0; i < persistenceDiagram1.size(); i++) {
        	vector<double> &point = persistenceDiagram1[i];
		costMatrix[i][persistenceDiagram2.size() + i] = 0;
		for(int p = 0; p < point.size(); p = p + 2){
			costMatrix[i][persistenceDiagram2.size() + i] =  max(costMatrix[i][persistenceDiagram2.size() + i], (abs(point[p+1]- point[p]))/2.0);
		}
		costMatrix[i][persistenceDiagram2.size() + i] = pow(costMatrix[i][persistenceDiagram2.size() + i], q);
	}
	i = 0;
        map<pair<double, double>, vector<vector<int> > >::iterator itr2;
  	for (int i = 0; i < persistenceDiagram2.size(); i++) {
        	vector<double> &point = persistenceDiagram2[i];
		costMatrix[persistenceDiagram1.size() + i][i] = 0;
		for(int p = 0; p < point.size(); p = p + 2){
			costMatrix[persistenceDiagram1.size() + i][i] =  max(costMatrix[persistenceDiagram1.size() + i][i], (abs(point[p+1]- point[p]))/2.0);
		}
		costMatrix[persistenceDiagram1.size() + i][i] = pow(costMatrix[persistenceDiagram1.size() + i][i], q);
	}

	for(int i = persistenceDiagram1.size(); i < numberOfPoints; i++){
		for(int j = persistenceDiagram2.size(); j < numberOfPoints; j++){
			costMatrix[i][j]= 0;//Matching diagonal points should take zero cost
		}
	}
	vector<int> Lmate;
	vector<int> Rmate;
	double cost = MinCostMatching(costMatrix, Lmate, Rmate);
	return cost;	
}

//Compute the Wasserstein distance between the persistence diagram "persistenceDiagram" and a persistence diagram consisting of infinitely many diagonal points (Each point in "persistence diagram" is matched to the closest diagonal point)
double vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::WassersteinDistance(vector<vector<double> > &persistenceDiagram){
	double distance = 0;
	for(int i = 0; i < persistenceDiagram.size(); i++){
		double cost=0;
		for(int p = 0; p < persistenceDiagram[i].size(); p = p + 2){
			cost = max(cost, abs(persistenceDiagram[i][p+1]-persistenceDiagram[i][p])/2.0);
		}
		distance += pow(cost, q); 
	}
	return distance;
}

//Compute the persistence of the point "point"
double vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::Persistence(vector<double> &point){
	double persistence = 0;
	for(int p = 0; p < point.size(); p = p + 2){
		persistence = max(persistence, abs(point[p + 1] - point[p]));
	}
	return persistence;
}


//Compute all possible Combinations of type strings
void vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::ConstructTypeStrings(string parentString, int cnt, vector<int> &types, vector<string> &typeStrings){
	if(cnt < this->Fields.size()){
		for(int t = 0; t < types.size(); t++){
			ConstructTypeStrings(parentString + DoubleToString(types[t]), cnt + 1, types, typeStrings);
		}
	}
	else{
		typeStrings.push_back(parentString);
	}
}


//Retrieve the string corresponding to a point "(a,b;c,d)" of the MDPD, representing the types (0-th ordinary/1st extended) of the persistence diagrams of the component Reeb graphs in the MDRG containing (a,b) and (c,d) 
string vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::GetTypeString(vector<int> types, vector<double> point){
	string typeString="";
	for(int i = 0; i < point.size(); i = i + 2){
		double birth = point[i];
		double death = point[i+1];
		if(birth <= death){
			typeString = typeString + DoubleToString(types[0]);
		}
		else{
			typeString = typeString + DoubleToString(types[1]);
		}
	}
	return typeString;
}

//Compute the Wasserstein Distance between the MDPDs corresponding to the JCNs "JCN1" and "JCN2"
double vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs::ComputeDistance(vtkGraph * JCN1, vtkGraph * JCN2){
	mdpdMap1 = *(new map<string, vector<vector<vector<double> > > >);
	mdpdMap2 = *(new map<string, vector<vector<vector<double> > > >);
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
 	JCN2->GetVertexData()->AddArray(nodeIdsJCN2);

	vtkMultiDimensionalReebGraph* mdrg1 = vtkMultiDimensionalReebGraph::New();
	vtkMultiDimensionalReebGraph* mdrg2 = vtkMultiDimensionalReebGraph::New();	
	for(int f = 0; f < this->Fields.size(); f++){
		mdrg1->AddField(this->Fields[f]);
		mdrg2->AddField(this->Fields[f]);
	}
	mdrg1->CreateMultiDimensionalReebGraph(JCN1);
	mdrg2->CreateMultiDimensionalReebGraph(JCN2);

	vector<vector<double> > combinedPersistenceDiagram1, combinedPersistenceDiagram2;
	vector<int> types;
	types.push_back(0);//0th-ordinary persistence diagram
	types.push_back(3);//1st-extended persistence diagram
	ConstructPersistenceDiagramMap(mdrg1, 0, 0, mdpdMap1, types, combinedPersistenceDiagram1, "0", 0);
	ConstructPersistenceDiagramMap(mdrg2, 0, 0, mdpdMap2, types, combinedPersistenceDiagram2, "0", 0);

	double wasserstein_distance_between_MDPDs = 0;
	std::set<string> fieldValues;
	for(map<string, vector<vector<vector<double> > > >::iterator it = mdpdMap1.begin(); it != mdpdMap1.end(); ++it){
	 	 fieldValues.insert(it->first);
	}
	for(map<string, vector<vector<vector<double> > > >::iterator it = mdpdMap2.begin(); it != mdpdMap2.end(); ++it){
	 	 fieldValues.insert(it->first);
	}
	string s="";
	vector<string> typeStrings;
	ConstructTypeStrings(s, 0, types, typeStrings);
	for(set<string>::iterator it = fieldValues.begin(); it != fieldValues.end(); ++it) {
		vector<vector<vector<double> > > *persistenceDiagramList1, *persistenceDiagramList2;
		 if(mdpdMap1.find(*it) != mdpdMap1.end()){
		 	persistenceDiagramList1 = &mdpdMap1[*it];
		 }
		else{
persistenceDiagramList1 = new vector<vector<vector<double> > >();
		}

		 if(mdpdMap2.find(*it) != mdpdMap2.end()){
		 	persistenceDiagramList2 = &mdpdMap2[*it];
		 }
		else{
persistenceDiagramList2 = new vector<vector<vector<double> > >();
		}

		vector<map<string,vector<vector<double> > > > persistenceDiagramList1Types((*persistenceDiagramList1).size());

		for(int i = 0; i < (*persistenceDiagramList1).size(); i++){
			vector<vector< double > > v;
			for(int t = 0; t < typeStrings.size(); t++){
				persistenceDiagramList1Types[i][typeStrings[t]]=v;
			}
			for(int j = 0; j < (*persistenceDiagramList1)[i].size(); j++){
				string typeString = GetTypeString(types, (*persistenceDiagramList1)[i][j]);
				persistenceDiagramList1Types[i][typeString].push_back((*persistenceDiagramList1)[i][j]);
			}					
		}

		vector<map<string,vector<vector<double> > > > persistenceDiagramList2Types((*persistenceDiagramList2).size());
		for(int i = 0; i < (*persistenceDiagramList2).size(); i++){
			vector<vector< double > > v;
			for(int t = 0; t < typeStrings.size(); t++){
				persistenceDiagramList2Types[i][typeStrings[t]]=v;
			}
			for(int j = 0; j < (*persistenceDiagramList2)[i].size(); j++){
				string typeString = GetTypeString(types, (*persistenceDiagramList2)[i][j]);
				persistenceDiagramList2Types[i][typeString].push_back((*persistenceDiagramList2)[i][j]);
			}
		}

		for(int t = 0; t < typeStrings.size(); t++){
			vector<vector<vector<double> > > persistenceDiagramsListByType1;
			vector<vector<vector<double> > > persistenceDiagramsListByType2;
			for(int i = 0; i < persistenceDiagramList1Types.size(); i++){
				persistenceDiagramsListByType1.push_back(persistenceDiagramList1Types[i][typeStrings[t]]);
			}
			for(int i = 0; i < persistenceDiagramList2Types.size(); i++){
				persistenceDiagramsListByType2.push_back(persistenceDiagramList2Types[i][typeStrings[t]]);
			}
			wasserstein_distance_between_MDPDs += WassersteinDistanceBetweenSetsOfPeristenceDiagrams(persistenceDiagramsListByType1, persistenceDiagramsListByType2);
		}
	}
	return pow(wasserstein_distance_between_MDPDs, 1.0/q);
}
