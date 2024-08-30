/*=========================================================================

 *	File: vtkDistanceBetweenMultiDimensionalReebGraphs.h
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/

// .NAME vtkDistanceBetweenMultiDimensionalReebGraphs - Computes the distance between the Multi-Dimensional Reeb Graphs (MDRGs) corresponding to the Joint Contour Nets (JCNs) of two multi-fields f and g 
//Input: Joint Countour Nets corresponding to bivariate fields f and g
//Output: Distance between the corresponding MDRGs

//Reference:
//[1] Yashwanth Ramamurthi and Amit Chattopadhyay. 2023. Topological Shape Matching using Multi-Dimensional Reeb Graphs. In Proceedings of the Thirteenth Indian Conference on Computer Vision, Graphics and Image Processing (ICVGIP '22). Association for Computing Machinery, New York, NY, USA, Article 5, 1–10.

#ifndef __vtkDistanceBetweenMultiDimensionalReebGraphs_h
#define __vtkDistanceBetweenMultiDimensionalReebGraphs_h
#include "vtkMutableDirectedGraph.h"
#include "vtkTree.h"
#include "vtkType.h"
#include<bits/stdc++.h>
#include "vtkMultiDimensionalReebGraph.h"
#include "vtkHandleDegeneracyInReebGraph.h"
#include "vtkPersistenceDiagramOfReebGraph.h"
using namespace std;

class VTK_META_EXPORT vtkDistanceBetweenMultiDimensionalReebGraphs : public vtkGraphAlgorithm
{
public:
        static vtkDistanceBetweenMultiDimensionalReebGraphs* New();
	vtkTypeMacro(vtkDistanceBetweenMultiDimensionalReebGraphs, vtkGraphAlgorithm);

	// Add the names of the component scalar fields of the bivariate fields (eg., f1,f2).
	void AddField(char *fieldName);

	//Compute the distance between the MDRGs corresponding to the JCNs "JCN1" and "JCN2" based on the persistence diagrams of a particular "type"
	//type 0 => 0-dimensional Persistence Diagram corresponding to the sublevelset filtration 
	//type 1 => 0-dimensional Persistence Diagram corresponding to the superlevelset filtration
	//type 2 => 0-dimensional Extended Persistence Diagram
	//type 3 => 1-dimensional Extended Persistence Diagram
 	//Note: Distances computed for types 0, 1, and 3 need to be multiplied by the weights w_0, w_1, and w_2, respectively, and the resulting values need to be be summed to obtain the distance d_T (see Equations 10 and 11 in [1])
	double ComputeDistance(vtkGraph * JCN1, vtkGraph * JCN2, int type);	
	
protected:
	//Segregation of the persistence diagrams of the component Reeb graphs of the MDRG of a multi-field f = (f1,f2,...,fn), into sets of the form S_{f_i}^{c}
	//persistenceDiagramMap[i][c] stores the list of persistence diagrams S_{f_i}^{c} (see Eqution (12) in [1])
	void ConstructPersistenceDiagramMap(vtkMultiDimensionalReebGraph *mdrg, int graphId, int fieldId, vector< map<double, vector<vector<pair<double,double> > > > > &persistenceDiagramMap, int type);

	//Computing the minimum-cost matching for the computation of the bottleneck distance (based on the Hungarian algorithm)
	double MinCostMatching(const vector<vector<double> > &cost, vector<int> &Lmate, vector<int> &Rmate);

	//Compute the bottleneck Distance between two persistence diagrams
	double BottleneckDistance(vector<pair<double, double> > persistenceDiagramOfReebGraph1, vector<pair<double, double> > persistenceDiagramOfReebGraph2);

	//Compute the bottleneck distance between the persistence diagram "persistenceDiagramOfReebGraph1" and a persistence diagram consisting of infinitely many diagonal points (Each point in "persistence diagram" is matched to the closest diagonal point)
	double BottleneckDistance(vector<pair<double, double> > persistenceDiagramOfReebGraph1);

	//Compute the sum of the distances between the sets of persistence diagrams in "persistenceDiagrams1" and "persistenceDiagrams2"
	//persistenceDiagrams1: consists of the persistence diagrams corresponding to the Reeb graphs in S_{f_{i+1}}^{c} (see Eqution (12) in [1])
	//persistenceDiagrams2: consists of the persistence diagrams corresponding to the Reeb graphs in S_{g_{i+1}}^{c} (see Eqution (12) in [1])
	double DistanceBetweenSetsOfPersistenceDiagrams(vector<vector<pair<double,double> > > &persistenceDiagrams1, vector<vector<pair<double,double> > > &persistenceDiagrams2);


	// Stores the names of the component scalar fields of the multi-fields (eg., f1,f2,..,fn).
	vector<char*> Fields;

	//Segregation of the persistence diagrams of the component Reeb graphs of the MDRGs based on the range values
	//persistenceDiagramMap1[i][c] stores the set of persistence diagrams S_{f_{i+1}}^{c} (see Eqution (12) in [1])
	//persistenceDiagramMap2[i][c] stores the set of persistence diagrams S_{g_{i+1}}^{c} (see Eqution (12) in [1])
	vector< map<double, vector<vector<pair<double,double> > > > > persistenceDiagramMap1, persistenceDiagramMap2;

    	vtkDistanceBetweenMultiDimensionalReebGraphs();
   	~vtkDistanceBetweenMultiDimensionalReebGraphs();

private:

};

#endif
