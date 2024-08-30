/*=========================================================================

 *	File: vtkHandleDegeneracyInReebGraph.h
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/

// .NAME vtkHandleDegeneracyInReebGraph - Handle degeneracy in the Reeb graph of a scalar field f 
//Input: A Reeb graph corresponding to a scalar field f
//Output: The degeneracy in the Reeb graph is removed


//Reference:
//[1] J. Tu, M. Hajij, and P. Rosen, “Propagate and pair: A single-pass approach to critical point pairing in reeb graphs,” in Advances in Visual Computing,
// Eds. Springer International Publishing, 2019.

#ifndef __vtkHandleDegeneracyInReebGraph_h
#define __vtkHandleDegeneracyInReebGraph_h
#include "vtkMultiDimensionalReebGraph.h" 
using namespace std;

class VTK_META_EXPORT vtkHandleDegeneracyInReebGraph : public vtkGraphAlgorithm
{
public:
        static vtkHandleDegeneracyInReebGraph* New();
 	vtkTypeMacro(vtkHandleDegeneracyInReebGraph, vtkGraphAlgorithm);

	/*Removes degeneracy in the Reeb graph by:
	1.  Handling degenerate critical nodes of the Reeb graph which violate the first Morse condition
	2.  Ensuring that no two critical nodes of the Reeb graph have the same critical value (second Morse condition)
	*/
	void HandleDegeneracy();

	//Sets the Reeb graph and the name of the associate scalar-field
	void SetParameters(vtkReebGraph* reebGraph, char *fieldName);

protected:

	//Breaks a non-degenerate minimum into a minimum and an up-fork
	void HandleDegenerateMinimum(int nodeID);
	
	//Breaks a non-degenerate maximum into a maximum and a down-fork
	void HandleDegenerateMaximum(int nodeID);
	
	//Break a degenerate up-fork with up-degree n into (n-1) up-forks
	void HandleDegenerateUpFork(int nodeID);

	//Break a degenerate down-fork with up-degree n into (n-1) down-forks
	void HandleDegenerateDownFork(int nodeID);
	
	//Break a double-fork into a down-fork and an up-fork
	void HandleDegenerateDoubleFork(int nodeID);
	
	//Ensures that no two critical nodes of the Reeb graph have the same critical value (second Morse condition)
	void AssignUniqueRangeValuesToCriticalNodes();
	
	vtkHandleDegeneracyInReebGraph();
	~vtkHandleDegeneracyInReebGraph();
	
	//Parameters for assigning range values to the nodes to satisfy the Morse conditions
	long double epsilon = 0.0001;
	int default_number_of_range_values = 10000;
	
	char *fieldName;
	vtkReebGraph *reebGraph;
};

#endif
