/*=========================================================================

 *	File: vtkJointContourNet.h
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/

// .NAME vtkJointContourNet - Compute the Joint Contour Net.
// .SECTION Description
// The Joint Contour Net (JCN) generalizes topological analysis of univariate scalar
// fields to multivariate fields.  The JCN is the corresponding multivariate extension 
// to the the Reeb Graph and Contour Tree: for univariate fields, the JCN yields a
// Contour Tree / Reeb Graph.  The JCN is constructed by quantizing the range space 
// of the input fields, yielding spatial regions (fragments) in which the range space
// values fall within a given interval.  Adjacent fragments that cover the same
// intervals for each variate are merged into slabs.  The JCN is then the dual graph
// over the slabs, where two slabs are connected if they share an edge in common.
//
// Two nodes of the JCN are connected iff:
// - all variables in the two corresponding slabs differ by at most 
//   one level of quantization;
// - at least one variable for which the slabs differs by one level 
//   of quantization;
//
// See:
// H. Carr, D. Duke, Joint Contour Nets: Multi-variate Analysis, Computation and Acceleration
// Proc. of IEEE Pacific Visualization 2013.
//
// The input to the VTK filter consists of an unstructured grid, the cells of which are 
// must [currently] be simplices.  There are three output ports. Port 0 is an undirected graph 
// defining the net structure.  Port 1 is geometry of fragments. Port 2 is geometry of slabs. 
// The multivariate fields over which the topology is defined must be 
// present as input point data: the vertex data of the output graph will contain a 
// copy of each such input field, with each array in the vertex data containing
// values quantized from the corresponding input field.  Output field data is always
// stored as double-precision values.
//
// .SECTION Caveats
// The filter assumes that cells are simplicial: each cell in in dataset whose domain 
// is R^M has M+1 points.
// Limits on the size of domain and range coordinates, and assumptions on numbers of
// quantization levels per dimension, are inhereted from vtkPolytopeGeometry.
//
// .SECTION See Also
// vtkPolytopeGeometry

#ifndef __vtkJointContourNet_h
#define __vtkJointContourNet_h

#include "vtkGraphAlgorithm.h"
#include "vtkMETAModule.h"
#include "vtkSmartPointer.h"
#include <vector>
#include "vtkPyramidTree.h"

class VTK_META_EXPORT vtkJointContourNet : public vtkGraphAlgorithm
{
public:
  static vtkJointContourNet *New();
  vtkTypeMacro(vtkJointContourNet,vtkGraphAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Flag to determine whether to use the input dataset's cell
  // structure (OFF), or a separate array on the input cell
  // attributes specifying the point IDs for each cell (ON).
  // Default is OFF.  If a separate array is used, the name of
  // the array should be set using the TopologyArrayName method.
  // NOTE: a separate array must currently be used for spatial 
  // domains of dimension > 3.
  vtkSetMacro(UseCellTopologyArray, bool);
  vtkGetMacro(UseCellTopologyArray, bool);
  vtkBooleanMacro(UseCellTopologyArray, bool);

  // Description:
  // Set/get the name of the array used for topological structure.
  vtkSetStringMacro(TopologyArrayName);
  vtkGetStringMacro(TopologyArrayName);

  // Description:
  // Flag to determine whether the domain (spatial) coordinates 
  // are carried by dataset points (OFF), or as a separate array 
  // in the input point data (ON).  If a separate coordinate array 
  // is // used, the array name must be specified using
  // SetSpatialCoordinateArrayName.  Default is OFF.
  // NOTE: as vtkPoints is limited to 3D, the separate array
  // must be used for spatial dimensions > 3.
  vtkSetMacro(UseSpatialCoordinateArray, bool);
  vtkGetMacro(UseSpatialCoordinateArray, bool);
  vtkBooleanMacro(UseSpatialCoordinateArray, bool);

  // Description:
  // Set/get the array name to use for spatial coordinates.
  vtkSetStringMacro(SpatialArrayName);
  vtkGetStringMacro(SpatialArrayName);

  // Description:
  // Determine whether to determine whether a slab lies on
  // the boundary (Default OFF).  If enabled, the filter
  // adds a vtkBitArray to the output vertex data,
  // where an entry is true iff the corresponding slab
  // has a facet that lies fully in the domain boundary.
  vtkSetMacro(IdentifyBoundarySlabs, bool);
  vtkGetMacro(IdentifyBoundarySlabs, bool);
  vtkBooleanMacro(IdentifyBoundarySlabs, bool);

  // Description:
  // The domain (spatial) dimension of the current dataset.
  vtkGetMacro(M, int);

  // Description:
  // The range (data) dimension of the current dataset.
  vtkGetMacro(N, int);
	
  // Description:
  // Add a field to be used in computing the JCN, and
  // specify the quantization level (width of slabs)
  // for this field.  The base of the field will be
  // computed from the input array range.
  void AddField(char *fieldName, double slabWidth);

  // Description:
  // Add a field to be used in computing the JCN, and
  // specify the quantization level (width of slabs)
  // for this field.  Also specify the base from which
  // slabs are computed.
  void AddField(char *fieldName, double slabWidth, double base);

  // Description:
  // Get the name of the ith field.
  char *GetField(int n);

  // Description:
  // Get the number of fields to be used for JCN construction.
  int GetNumberOfFields();

  // Description:
  // Clear all stored fields.
  void ClearFields();

  // Description:
  // Get the JCN graph.
  vtkUndirectedGraph *GetTopologyGraph();
	
  // Description:
  // Get the ith output.  Currently only one output is defined,
  // but class may be extended to generate slab/fragment geometry.
  vtkDataObject *GetOutput(int index);

  // Description:
  // Determine whether quantization starts from the minimum
  // field value, or by centering the first "slab" over the
  // mimumum field value. Default is OFF.
  vtkSetMacro(CenterFirstSlab, bool);
  vtkGetMacro(CenterFirstSlab, bool);
  vtkBooleanMacro(CenterFirstSlab, bool);
	
protected:
  vtkJointContourNet();
  ~vtkJointContourNet();

  // Override to specify support for any vtkDataSet input type.
  virtual int FillInputPortInformation(int port, vtkInformation* info);
  virtual int FillOutputPortInformation(int port, vtkInformation* info);
	
  virtual int RequestDataObject(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
	
  // Main implementation.
  virtual int RequestData(vtkInformation*,
                          vtkInformationVector**,
                          vtkInformationVector*);

	// Union-find structure used to compute fragment equivalence.
  vtkSmartPointer<vtkIdTypeArray> UF;		
  vtkIdType Find(vtkIdType x);

  // Global parameters.
  bool CenterFirstSlab;
  bool UseCellTopologyArray;
  bool UseSpatialCoordinateArray;
  bool IdentifyBoundarySlabs;

  char *TopologyArrayName;
  char *SpatialArrayName;

  int M;  // domain (spatial) dimension
  int N;  // range (data) dimension

  // Set of fields that define the multifield structure.		
  std::vector<struct _Field*> Fields;

  // Utility function for converting a value v in field f into
  // the corresponding slab coordinate given the quantization
  // settings.
  double SlabCoordinate(int f, double v);

private:
  vtkJointContourNet(const vtkJointContourNet&);    // Not implemented.
  void operator=(const vtkJointContourNet&);  // Not implemented.
};

#endif

