/*=========================================================================

 *	File: vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs.h
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/

// .NAME vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs - Compute the Wasserstein Distance between the Multi-Dimensional Persistence persistence diagrams (MDPDs) of bivariate fields f and g 
//Input: Joint Countour Nets corresponding to bivariate fields f and g
//Output: Wasserstein Distance between the MDPDs

//Reference:
//[1] Yashwanth Ramamurthi and Amit Chattopadhyay, "A Topological Distance Between Multi-Fields Based on Multi-Dimensional Persistence Diagrams," in IEEE Transactions on Visualization and Computer
//Graphics, vol. 30, no. 9, pp. 5939-5952, Sept. 2024


#ifndef __vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs_h
#define __vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs_h
#include "vtkMutableDirectedGraph.h"
#include "vtkTree.h"
#include "vtkType.h"
#include<bits/stdc++.h>
#include "vtkMultiDimensionalReebGraph.h"
#include "vtkHandleDegeneracyInReebGraph.h"
#include "vtkPersistenceDiagramOfReebGraph.h"
using namespace std;

class VTK_META_EXPORT vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs : public vtkGraphAlgorithm
{
public:
        static vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs* New();
	vtkTypeMacro(vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs, vtkGraphAlgorithm);

	// Add the names of the component scalar fields of the bivariate fields (eg., f1,f2).
	void AddField(char *fieldName);

	// Set the parameter q for computing the Wasserstein distance (see Equations 10 and 11 in [1])
	void SetParameterQ(double q);

	//Compute the Wasserstein Distance between the MDPDs corresponding to the JCNs "JCN1" and "JCN2"
	double ComputeDistance(vtkGraph * JCN1, vtkGraph * JCN2);

	//Convert a double (floating point) value to a string
	string DoubleToString(double number);
	
	//The parameter q in the equation for Wasserstein distance (see Equations 10 and 11 in [1])
	double q = 2.0;
	
protected:
	
	//Computes the Multi-Dimensional Persistence Diagram PD_{f} corresponding to an MDRG MR_{f}. 
	// The points in the MDPD are stored in groups of the form PD_p^{x_1}(f)//
	// For a real-number v, multiDimensionalPersistenceDiagramMap[v] stores a list of {PD_p^{x_1}(f): \overline{f_1}(p) = val} (see Equation 4 in [1] for the definition of PD_p^{x_1}(f))
	// The segregation of the points in the MDPD into groups helps in the checking of criteria (C1)-(C3) while computing the Wasserstein distance between two MDPDs
	void ConstructPersistenceDiagramMap(vtkMultiDimensionalReebGraph *mdrg, int graphId, int fieldId, map<string, vector<vector<vector<double> > > > &multiDimensionalPersistenceDiagramMap, 	vector<int> &types, vector<vector<double> > combinedParentPersistenceDiagram, string fieldValue, double parentFieldValue);

	//Computing the minimum-cost matching for the computation of the Wasserstein distance (based on the Hungarian algorithm)
	double MinCostMatching(const vector<vector<double> > &cost, vector<int> &Lmate, vector<int> &Rmate);

	//Compute the Wasserstein distance between sets of points in the MDPD grouped into different categories based on the matching criteria (C1)-(C3)
	double WassersteinDistanceBetweenSetsOfPeristenceDiagrams(vector<vector<vector<double> > > &persistenceDiagrams1, vector<vector<vector<double> > > &persistenceDiagrams2);

	//Compute the Wasserstein Distance between two persistence diagrams
	double WassersteinDistance(vector<vector<double> > &persistenceDiagram1, vector<vector<double> > &persistenceDiagram2);

	//Compute the Wasserstein distance between the persistence diagram "persistenceDiagram" and a persistence diagram consisting of infinitely many diagonal points (Each point in "persistence diagram" is matched to the closest diagonal point)
	double WassersteinDistance(vector<vector<double> > &persistenceDiagram);

	//Compute the persistence of the point "point"
	double Persistence(vector<double> &point);

	//Compute all possible Combinations of type strings
	void ConstructTypeStrings(string parentString, int cnt, vector<int> &types, vector<string> &typeStrings);

	//Retrieve the string corresponding to a point "(a,b;c,d)" of the MDPD, representing the types (0-th ordinary/1st extended) of the persistence diagrams of the component Reeb graphs in the MDRG containing (a,b) and (c,d)
	string GetTypeString(vector<int> types, vector<double> point);

	// Stores the names of the component scalar fields of the bivariate fields (eg., f1,f2).
	vector<char*> Fields;

	//Stores the two MDPDs
	map<string, vector<vector<vector<double> > > > mdpdMap1, mdpdMap2;

    	vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs();
   	~vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs();

private:

};

#endif
