/*=========================================================================

 *	File: vtkPolytopeGeometry.h
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/

// .NAME vtkPolytopeGeometry - store and process polytopes in R^m->R^n.
// 
// .SECTION Description
// vtkPolytopeGeometry provides methods for creating, manipulating and
// saving (in a compact form) polytopes whose vertices are points in 
// an arbitrary dimensional domain space (R^m) and where each vertex is
// associated with an n-dimensional sample (in R^n).
//
// Two levels of storage are utilized.  Polytopes are initially created
// using an explicit representation:
// 1.  coordinates in R^m are stored uniquely in a coordinate table;
// 2.  each facet of a polytope (edge, polygon, etc) is explicitly
//     stored, and assigned a unique index.  Facets common to more than
//     one sub-polytope are thus shared.
// 3.  each range sample is stored, although not uniquely, i.e. if two
//     vertices have the same range sample values, those two sets of 
//     values will be stored distinctly.
//
// Operations are provided to:
// * add a coordinate to the active store
// * add a polytope (which must be built from existing polytopes and/or)
//   coordinates in the active store.
// * split a polytope into a collection of fragments by slicing the
//   polytope against a set of cutting hyper-planes orthogonal to a 
//   given axis.
// * retrieve the minimum (wrt bounding box) coordinate of a polytope.
// * copy active polytopes into storage.
//
// Results of fragmenting polytopes are cached, as sub-polytopes will 
// typically be shared between polytope facets and/or multiple polytopes.
//
// NOTE: although vertices are technically 0-topes, for efficiency 
// vertices are not stored explicitly as polytopes, only as indices
// into a vertex table.  Although the interpretation of an id as a
// polytope or vertex index is usually clear from context, we also
// separate the two sets by mapping vertex indices to negative values,
// using the formula vertex-index = -1 - vertex-offset.
//
// Inactive (stored) polytopes are currently represented by 
// - the minimum coordinate of the polytope's bounding box
// - the center point of each facet. 
// - the ids of the polytopes adjacent to each facet center 
//   point.  Where a facet lies on a boundary, the first of
//   the ids is the facet's polytope, the second is "-1".
//
// Once stored, a polytope cannot be made active.
//
// .SECTION Caveats
// Design and implementation of this class has been driven by the
// implementation of vtkJointContourNet.  There is considerable
// scope for generalization, which may result in interface changes.
//
// Storage of non-active polytopes assumes that at most two 
// m-dimensional polytopes share an (m-1)-dimensional facet.
//
// Program limits (see constants below) are not checked at runtime.


#ifndef __vtkPolytopeGeometry_h
#define __vtkPolytopeGeometry_h

#include "vtkIdTypeArray.h"
#include "vtkMETAModule.h"
#include "vtkObject.h"
#include "vtkPoints.h"
#include "vtkPyramidTree.h"
#include "vtkSmartPointer.h"
#include <vector>


// Maximum number of components for domain and range coordinates.
const unsigned int MaxRangeDim = 20;
const unsigned int MaxDomainDim = 20;

// Maximum number of slices per polytope.
const unsigned int MaxNrSlicesPerPtope = 256;

// Bounds on polytope structure.
const int MaxNrFacetsPerPtope = 50;
const int MaxNrPointsPerPtope = 5000;

class vtkMergePoints;
class vtkDoubleArray;
class vtkPoints;

// BTX
// Record structure for active polytopes.
struct _Ptope {
    vtkIdType start;      // offset of first facet id in Facet array.
    short size;           // number of facets
    short dim;            // spatial dimensionality
    double *center;       // if non-null, the center of the polytope.
    double *minCoord;     // if non-null, the lowest corner of bounding box
    int splitDim;         // dimension on which any cached fragments were generated
    vtkIdType *fragments; // if non-null, cached fragments
    vtkIdType *ends;      // if non-null, cached ends from cutting planes
    };
// ETX

class VTK_META_EXPORT vtkPolytopeGeometry : public vtkObject {
public:
  // Description:
  // Standard VTK interface for creating and printing an object.
  static vtkPolytopeGeometry *New();
  vtkTypeMacro(vtkPolytopeGeometry, vtkObject);
  virtual void PrintSelf(ostream &os, vtkIndent indent);

  // Description:
  // Initialize a geometry object, specifying domain and range dimensions,
  // and providing in advance.
  void Initialize(const short m, const short n, const double *domainBounds);  

  // Description:
  // Clear the active polytope store to empty, and initialise it to
  // receive new polytopes within the specified domain (spatial)
  // bounds.  "domainBounds" points to 2M doubles containing the
  // min & max bounds for each of the M ordinates.
  void ResetForNextCell(const double *domainBounds);

  // Description:
  // Specialized version of PrintSelf, which prints the top
  // level (polytope id and list of points, domain and range
  // coordinates) for the "sz" active polytopes whose indices
  // are pointed to bt "top".
  void PrintToplevel(ostream &os, vtkIdType *top, int sz);

  // Description:
  // Copy information on active polytopes into the store.
  // Active polytopes are not changed/deleted.
  void StoreActivePolytopes(vtkIdType *ids, int nr);

  // Description:
  // Add a simplex specified by a list of point ids, and
  // a dimension, to the active polytopes.  The point ids
  // must refer to coordinates that have already been added.
  vtkIdType BuildFromSimplex(vtkIdType *pt, int d);

  // Description:
  // Add a vertex (0-d polytope), specified by a coordinate in 
  // domain space, and a coordinate in range space.
  vtkIdType AddCoord(const double *dc,const double *rc);

  // Description:
  // Add a polytope to the active polytope set.  The polytope
  // is specified by a dimensionality, and a list of "sz" ids
  // to the polytope/vertex ids of its facets.
  vtkIdType AddPtope(vtkIdType *pt, int d, int sz);

  // Description:
  // Access the domain coordinate of an active polytope,
  // returning coordinate in buffer provided by caller.
  void GetActiveDomCoord(vtkIdType id, double *c);

  // Description:
  // Access the range coordinate of an active polytope,
  // returning coordinate in buffer provided by caller.
  void GetActiveRngCoord(vtkIdType id, double *c);

  // Description:
  // Return a specific component of the domain coordinate
  // of an active polytope.
  double GetActiveDomCoordComponent(vtkIdType id, int c);

  // Description:
  // Return a specific component of the range coordinate
  // of an active polytope.
  double GetActiveRngCoordComponent(vtkIdType id, int c);

  // Description:
  // Return a specific component of the range coordinate
  // of a stored polytope.
  double GetStoredRangeComponent(vtkIdType id, int c);

  // Description:
  // Split the active polytope with id "ptid" into a set 
  // of polytope fragments by slicing ptid against "nrts" 
  // threshold values.  Each threshold value specifies a
  // cutting plane orthogonal to axis "field"; thresholds
  // values must appear in increasing value.
  // Each intersections of a cutting plane and the polytope
  // also result in a lower-dimensional 'end' facet.  The
  // fragments and ends are added to the active polytopes, 
  // and their ids returned in the caller-allocated 
  // "frags" and "ends" buffers.  
  // *NOTE*: the buffers must
  // contain at least one more slot than the number of
  // cutting planes: the ith slot in frags holds the 
  // fragment that is LESS THAN the ith cutting plane;
  // the "nrts" slot holds the residue of the polytope 
  // that is greater than the highest cutting plane.
  void Split(
    vtkIdType ptid, 
    double *thresholds, 
    int nrts, 
    int field, 
    vtkIdType *frags, 
    vtkIdType *ends);

  // Description:
  // Retrieve the object holding the center point of
  // each stored polytope facet.
  //vtkGetObjectMacro(CenterPoints, vtkPoints);
  vtkGetObjectMacro(CenterLocator,vtkPyramidTree);

  // Description:
  // Retrieve the two-component array giving the ids
  // of the one or two polytopes that share a facet.  
  vtkGetObjectMacro(CenterFacets, vtkIdTypeArray);
  vtkGetObjectMacro(CenterFacetsId, vtkIdTypeArray);

  // Description:
  // Return the number of stored polytopes.
  vtkGetMacro(NrStoredPtopes, vtkIdType);
  
  // Descriptioin:
  // Return the facet starting index and size, given a facet id. 
  void GetPtopeFacets(vtkIdType id, vtkIdType *&ids, int &nr);
  
protected:
  vtkPolytopeGeometry();
  ~vtkPolytopeGeometry();
  
private:
  vtkPolytopeGeometry(const vtkPolytopeGeometry&); // Not implemented
  void operator=(const vtkPolytopeGeometry&); // Not implemented

  static const int TABSIZE  =      9001;
  static const int SETSIZE  =        32;
  static const int MaxNrPtopes =    9999;
  static const int MaxNrFacets =  50000;
  
  // dimensionality
  short M;    // domain (spatial) dimension
  short N;    // range (data) dimension

  // Active polytope storage.  
  // Polytopes are held in an array of records indexed directly by
  // polytope ID.  For algorithms that require fast lookup of polytopes
  // based on structure, a hashtable is used, hashing on the *SET* of
  // points defining the polytope.  The hash key is computed from the
  // *sorted* list of polytope components (points or facet polytope ids).
  // For efficiency (avoiding repeated memory allocation) the facet ids
  // of all polytopes are held in a single extensible array, with each
  // polytope holding a start offset and number of facets.
  struct _Ptope Ptope[MaxNrPtopes];
  vtkIdType Facet[MaxNrFacets];

  // Hash table:
  // We use open hashing, implemented by three structures:
  // - Index records the polytope id for a given hash.
  // - SetSize records the key size (i.e. number of polytope components)
  //   for occupied entries, and -1 for unused entries.  It thuse doubles
  //   as an occupancy map and fast reject test for matches against a key.
  // - PointSet holds the key for each hash entry.  
  //   CURRENTLY ASSUMES FIXED UPPER-BOUND (SETSIZE) ON KEY SIZE.

  vtkIdType Index[TABSIZE];
  short SetSize[TABSIZE];
  vtkIdType PointSet[TABSIZE][SETSIZE];

  // Global counters for active polytope and facet tables.
  vtkIdType NextFacet;
  vtkIdType NextPtope;

  // index for stored polytopes
  vtkIdType NrStoredPtopes;

  // Buffer for recording polytope centers.
  double Center[MaxDomainDim];

  // Helper functions for Split method, handling
  // specific cases.
  void SplitEdge(
    struct _Ptope *entry, 
    double *thresholds, 
    int nrts, 
    int field,
    vtkIdType *frags, 
    vtkIdType *ends);
    
  void SplitPtope(
    struct _Ptope *entry, 
    double *thresholds, 
    int nrts, 
    int field,
    vtkIdType *frags, 
    vtkIdType *ends);

  // Helper function for Split method, extracting a polytope
  // of a specific dimensionality from candidate set either by
  // finding a single instance of the required size, or by
  // assembling a candidate from a set of lower-dimensional facets.
  vtkIdType Select(
    vtkIdType *pts, 
    int nrts,
    int nrps,
    int slice,
    int startDim, int endDim);

  double *GetPtopeCoordMin(vtkIdType id);
  //void GetPtopeFacets(vtkIdType id, vtkIdType *&ids, int &nr);
  double *GetPtopeCenter(vtkIdType id);
  
  // Acceleraton strutures for active polytopes
  vtkSmartPointer<vtkPyramidTree> ActiveDomCoord;
     
  // Structures for stored polytopes.
  vtkSmartPointer<vtkPyramidTree> StoredDomCoord;
  

  // Storage and accelerated lookup for facet center points.
  vtkPyramidTree *CenterLocator;


  // Map linking polytopes to facet centers.
  // Each m-dimensional polytope has a set of 
  // (m-1)-dimensional facets.  We assume each facet 
  // is shared by at most two polytopes.  The two 
  // components in the ith entry of CenterFacets are 
  // the point ids of the one (or two) polytopes 
  // containing the ith in the CenterPoints array
  // as defined above.
  vtkIdTypeArray *CenterFacets;
  vtkIdTypeArray *CenterFacetsId;

  // Storage for range coordinates (multi-field scalars).
  vtkSmartPointer<vtkDoubleArray> ActiveScalars;
  vtkSmartPointer<vtkDoubleArray> StoredScalars;
};

#endif
