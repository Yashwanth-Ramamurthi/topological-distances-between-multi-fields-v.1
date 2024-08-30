/*=========================================================================

 *	File: vtkPersistenceDiagramOfReebGraph.h
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/

// .NAME vtkPersistenceDiagramOfReebGraph - Compute the persistence diagrams corresponding to the Reeb graph $\RG_{f}$ of a scalar field f. 
//Reference:
//Ulrich Bauer, Xiaoyin Ge, and Yusu Wang. 2014. Measuring Distance between Reeb Graphs. In Proceedings of the thirtieth annual symposium on Computational geometry (SOCG'14).
//Association for Computing Machinery, New York, NY, USA

//See also the following article for the construction of different persistence diagrams:
//Yashwanth Ramamurthi and Amit Chattopadhyay, "A Topological Distance Between Multi-Fields Based on Multi-Dimensional Persistence Diagrams," in IEEE Transactions on Visualization and Computer Graphics, vol. 30, no. 9, pp. 5939-5952, Sept. 2024

#ifndef __vtkPersistenceDiagramOfReebGraph_h
#define __vtkPersistenceDiagramOfReebGraph_h
#include "vtkMutableDirectedGraph.h"
#include "vtkType.h"
#include<bits/stdc++.h>
#include "vtkMultiDimensionalReebGraph.h"

using namespace std;

class VTK_META_EXPORT vtkPersistenceDiagramOfReebGraph : public vtkGraphAlgorithm
{
public:
        static vtkPersistenceDiagramOfReebGraph* New();
	vtkTypeMacro(vtkPersistenceDiagramOfReebGraph, vtkGraphAlgorithm);

	//Obtain the persistence diagram of a specific type
	//type 0 => 0-dimensional Persistence Diagram corresponding to the sublevelset filtration 
	//type 1 => 0-dimensional Persistence Diagram corresponding to the superlevelset filtration
	//type 2 => 0-dimensional Extended Persistence Diagram
	//type 3 => 1-dimensional Extended Persistence Diagram
     	vector<pair<double, double> > GetPersistenceDiagram(int type);


	//Constructs the four types of persistence diagrams
     	void ConstructPersistenceDiagram();

	vtkReebGraph* reebGraph;
	char *fieldName;	
protected:
	//Returns the id of the node with the minimum range (f) value by traversing the Reeb graph from nodeId to all the nodes with range (f) value  at most "maximumRangeValue"
	vtkIdType GetMinimumNode(int nodeId, double maximumRangeValue, vector<int> &path);

	//Returns the id of the node with the maximum range (f) value by traversing the Reeb graph from nodeId to all the nodes with range (f) value  at least "minimumRangeValue"
	vtkIdType GetMaximumNode(int nodeId, double minimumRangeValue, vector<int> &path);

	//Returns the id of the node with the maximum range (f) value by traversing the Reeb graph from nodeId to all the nodes with range (f) value  at least "minimumRangeValue"
	void VisitNodesWithinASpecifiedRange(int nodeId, double min_range, double max_range, vector<bool> &visited);
	
	//Stores the persistence diagram corresponding to the sublevelset filtration (with the infinity point (a,\infty) replaced with (a, global maximum of f))
	vector<pair<double, double> > persistenceDiagramSublevelSetFiltration;
	
	//Stores the persistence diagram corresponding to the superlevelset filtration
	vector<pair<double, double> > persistenceDiagramSuperLevelSetFiltration;
	
	//Stores the 0th-extended persistence diagram
	vector<pair<double, double> > extendedPersistenceDiagram_dimension_0;	
	
	//Stores the 1st-extended persistence diagram
	vector<pair<double, double> > extendedPersistenceDiagram_dimension_1;	

    	vtkPersistenceDiagramOfReebGraph();
   	~vtkPersistenceDiagramOfReebGraph();

private:

};

#endif
