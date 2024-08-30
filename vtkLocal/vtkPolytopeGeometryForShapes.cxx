/*=========================================================================

 *	File: vtkPolytopeGeometryForShapes.cxx
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
#include "vtkDoubleArray.h"
#include "vtkIndent.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkMergePoints.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPolytopeGeometryForShapes.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include <vector>
#include <queue>

// Tolerance used when comparing components of points.
#define EPSILON 1.0e-10

vtkStandardNewMacro(vtkPolytopeGeometryForShapes);

//----------------------------------------------------------------------------
// Polytope storage 
//----------------------------------------------------------------------------

vtkPolytopeGeometryForShapes::vtkPolytopeGeometryForShapes()
{
  // Acceleration structures for active polytopes
  this->ActiveDomCoord     = vtkSmartPointer<vtkPyramidTree>::New();
    
  // Initialize active polytope table.
  // Ptope entry 0 is reserved as a "null" (-1-dimensional) ptope.
  this->Ptope[0].dim   = -1;
  this->Ptope[0].start =  0;
  this->Ptope[0].size  =  0;
  this->Ptope[0].minCoord = NULL;
  this->Ptope[0].fragments = NULL;
  this->Ptope[0].ends = NULL;
  this->NextPtope = 1;

  // Counter for stored (inactive) polytopes.  
  this->NrStoredPtopes = 0;

  this->CenterFacets = vtkIdTypeArray::New();
  this->CenterFacets->Register(this);
  this->CenterFacets->SetNumberOfComponents(2);  
  
  this->CenterFacetsId = vtkIdTypeArray::New();
  this->CenterFacetsId->Register(this);
 
  // Stored point table and lookup
  this->StoredDomCoord= vtkSmartPointer<vtkPyramidTree>::New();

  // Storage for facet center, and lookup.
  this->CenterLocator = vtkPyramidTree::New();
  this->CenterLocator->Register(this);
}

vtkPolytopeGeometryForShapes::~vtkPolytopeGeometryForShapes() 
{
  this->CenterLocator->Delete();
  this->CenterFacets->Delete();
  this->CenterFacetsId->Delete();
}

// Initialize structures for given domain (M) and range (N) 
// dimensionality.
void vtkPolytopeGeometryForShapes::Initialize(
  const short m, 
  const short n,
  const double *domainBounds)
{

  // Initialize maximum dimensions of polytopes & data to be stored.
  this->M = m;
  this->N = n;
  
  if (this->M == 2)
    {
    this->CenterFacetsId->SetNumberOfComponents(4);  
    }
  else 
    {
    this->CenterFacetsId->SetNumberOfComponents(2);  
    }
    
  this->ActiveDomCoord->Initialize(this->M,domainBounds);

  // Initialize structure for accelerating looking up of 
  // active domain points.;
  this->StoredDomCoord->Initialize(this->M,domainBounds);

  // Initialize structure for accelerating looking up of 
  // active facet center points.
  this->CenterLocator->Initialize(this->M,domainBounds);
  

  // Create data arrays to hold range values for both active 
  // and stored polytopes.
  this->ActiveScalars = vtkSmartPointer<vtkDoubleArray>::New();
  this->StoredScalars = vtkSmartPointer<vtkDoubleArray>::New();
  this->ActiveScalars->SetNumberOfComponents(this->N);
  this->StoredScalars->SetNumberOfComponents(this->N);
    
  for (int i = 1; i < MaxNrPtopes; i++)
     {
     this->Ptope[i].center = NULL;
     this->Ptope[i].minCoord = NULL;
     this->Ptope[i].fragments = NULL;
     this->Ptope[i].ends = NULL;
     }
}

// ------------------------------------------------------------------

// Clear any active polytopes, resetting active storage structures,
// and re-initialize the structures for new polytopes within the
// specified spatial (domain) bounds.
void vtkPolytopeGeometryForShapes::ResetForNextCell(const double *domainBounds)
{
  struct _Ptope *entry = this->Ptope+1;
  
  // Clear polytope table.
  // Entry 0, the Nil polytope, is unaffected.
  for (int i = 1; i < this->NextPtope; i++, entry++)
    {    
    if (entry->center)    delete [] entry->center;
    if (entry->minCoord)  delete [] entry->minCoord;
    if (entry->fragments) delete [] entry->fragments; 
    if (entry->ends)      delete [] entry->ends; 
    }

    
  // Reset polytope and facet tables counters.
  this->NextFacet = 0;
  this->NextPtope = 1;
  
  // Reset hash-table support.
  for (int i = 1; i < vtkPolytopeGeometryForShapes::TABSIZE; i++)
    {
    this->SetSize[i] = -1;
    }
  
  // Reset domain coordinate storage and lookup acceleration.
  this->ActiveDomCoord->Initialize(this->M,domainBounds);
    
}

// ------------------------------------------------------------------

// Add a polytope specified by a list of 'sz' point/sub-tope indices.

vtkIdType vtkPolytopeGeometryForShapes::AddPtope(vtkIdType *pt, int d, int sz)
{
  vtkIdType p;
  vtkIdType temp[MaxNrPointsPerPtope];
  unsigned long long h;
  struct _Ptope *entry;
  double fp, fq;
  double c1[MaxRangeDim], c2[MaxRangeDim];

  vtkIdType a, b, c; 

  // Copy component ids to avoid interference with caller's data, and
  // sort ids to provide a unique hash key for a given set of ids.
  for (int i = 0; i < sz; i++) 
    {
    temp[i] = pt[i];
    }
  std::sort(temp, temp+sz);

  // Hash the key.  Hash function to use is determined by number of indices
  // in the polytope, to try and get a wider spread of values over the common
  // low-dimensional cases.
   switch (sz) {
    case 1: 
      h = temp[0]; 
      break;
    case 2:
      h = (temp[1] << 16) + temp[0];
      break;
    case 3:
      h = (temp[2] << 24) + (temp[1] << 16) + temp[0];
      break;
    default:
      h = 0;
      for (int j = 0; j < sz; j++) { h += temp[j]; h = h << 4; }
    }

  // compute index into hash table.
  
  h = h % vtkPolytopeGeometryForShapes::TABSIZE;

  // linear probing: starting from initial hash location, probe
  // consecutive locations until either an empty slot is found,
  // or an exact match is found.  In the former case, create 
  // and initialize a new polytope.

  while (true)
    {
    if (this->SetSize[h] < 0)
      {
      // Empty slot, poltope not present.
      // Copy polytope key into this slot.
      this->SetSize[h] = sz;
      for (int i = 0; i < sz; i++)
        {
        this->PointSet[h][i] = temp[i];
        }
      // Initialize the polytope record.
      entry = this->Ptope + this->NextPtope;
      entry->start = this->NextFacet;
      entry->dim   = d;
      entry->size  = sz;
      entry->center = NULL;
      entry->fragments = NULL;
      entry->ends = NULL;
      // Copy facet indices into facet array,
      for (int i = 0; i < sz; i++)
        {
        this->Facet[this->NextFacet++] = temp[i];
        }
      
       entry->minCoord = new double[this->N];
      
      // Compute minimum coordinate of polytope.
      if (d == 1)
        {
        // Edge: minimum is componentwise min over 
        // coordinates of end points.
        this->GetActiveRngCoord(pt[0], c1);
        this->GetActiveRngCoord(pt[1], c2);
        for (int j = 0; j < this->N; j++)
          {
          entry->minCoord[j] = std::min(c1[j], c2[j]);
          }
        }
      else
        {
        // Non-edge: Initialise minimum to min coordinate of
        // first component, then update against min of each
        // subsequent component.
        for (int j = 0; j < this->N; j++)
          {
          entry->minCoord[j] = this->Ptope[pt[0]].minCoord[j];
          }
        for (int i = 1; i < sz; i++)
          {
          for (int j = 0; j < this->N; j++)
            {
            entry->minCoord[j] = std::min( 
                entry->minCoord[j],
                this->Ptope[pt[i]].minCoord[j] );
            }
          }
        }
      
      // Update hash index entry to point to this polytope.
      this->Index[h] = this->NextPtope++;
      break;
      }
    else 
      {
      // Hash table contains an entry at this position.
      // Quick test: does the size match?
      if (this->SetSize[h] == sz)
        {
        // Test whether each component matches.
        bool equal = true;
        for (int j = 0; equal && j < sz; j++)
          {
          equal &= this->PointSet[h][j] == temp[j];
          }
        // If matched, the polytope is present at index h.
        if (equal) 
          {
          break;
          }
        } 
      // No match: advance to next slot in table.
      h++;
      if (h == vtkPolytopeGeometryForShapes::TABSIZE) h = 0;
      }
    }

  // Guaranteed to have found or created the polytope in slot h.
  return this->Index[h];
} // AddPtope


// ------------------------------------------------------------------

// Construct an active polytope from a set of coordinates defining
// a d-dimensional simplex.

vtkIdType vtkPolytopeGeometryForShapes::BuildFromSimplex(vtkIdType *pt, int d)
{
  vtkIdType points[MaxNrFacetsPerPtope];
  vtkIdType facets[MaxNrFacetsPerPtope];
  vtkIdType lines[3];
  vtkIdType triangles[4];
  
  int size = d+1;
  vtkIdType sid;

  switch(d) {
    
    case 0:     // must be a coordinate, already in table.
      break;
    
    case 1:     // edge defined by two points, use directly.
      sid = this->AddPtope(pt, d, d+1);
      break;

    case 2:     // triangle: build three edge polytopes, then define a polygon.
      points[0] = pt[0]; points[1] = pt[1]; lines[0] = this->AddPtope(points, 1, 2);
      points[0] = pt[0]; points[1] = pt[2]; lines[1] = this->AddPtope(points, 1, 2);
      points[0] = pt[1]; points[1] = pt[2]; lines[2] = this->AddPtope(points, 1, 2);
      sid = this->AddPtope(lines, 2, 3);
      break;
      
    case 3:     // tetrahedron: build four triangles, one per face.
      // Face 0 1 2
      points[0] = pt[0]; points[1] = pt[1]; lines[0] = this->AddPtope(points, 1, 2);
      points[0] = pt[0]; points[1] = pt[2]; lines[1] = this->AddPtope(points, 1, 2);
      points[0] = pt[1]; points[1] = pt[2]; lines[2] = this->AddPtope(points, 1, 2);
      triangles[0] = this->AddPtope(lines, 2, 3);
    
      // Face 0 1 3
      points[0] = pt[0]; points[1] = pt[1]; lines[0] = this->AddPtope(points, 1, 2);
      points[0] = pt[0]; points[1] = pt[3]; lines[1] = this->AddPtope(points, 1, 2);
      points[0] = pt[1]; points[1] = pt[3]; lines[2] = this->AddPtope(points, 1, 2);
      triangles[1] = this->AddPtope(lines, 2, 3);

      // Face 0 2 3
      points[0] = pt[0]; points[1] = pt[2]; lines[0] = this->AddPtope(points, 1, 2);
      points[0] = pt[0]; points[1] = pt[3]; lines[1] = this->AddPtope(points, 1, 2);
      points[0] = pt[2]; points[1] = pt[3]; lines[2] = this->AddPtope(points, 1, 2);
      triangles[2] = this->AddPtope(lines, 2, 3);

      // Face 1 2 3
      points[0] = pt[1]; points[1] = pt[2]; lines[0] = this->AddPtope(points, 1, 2);
      points[0] = pt[1]; points[1] = pt[3]; lines[1] = this->AddPtope(points, 1, 2);
      points[0] = pt[2]; points[1] = pt[3]; lines[2] = this->AddPtope(points, 1, 2);
      triangles[3] = this->AddPtope(lines, 2, 3);
      sid = this->AddPtope(triangles, 3, 4);
      break;
      
    default:
      // Reduce dimensionality by one, and recurse.
      // Given n points, form n (d-1)-dimensional facets by
      // systematically omitting one point at a time.
      for (int omit = 0; omit < size; omit++)
        {
        for (int i = 0, p = 0; i < size; i++)
          {
          if (i == omit) continue;
          points[p++] = pt[i];
          }
        facets[omit] = this->BuildFromSimplex(points, d-1);
        } // for each face
      sid = this->AddPtope(facets, d, d+1);        

     } // switch dimension

  // return the id of the top-level polytope.
  return sid;
}

// ------------------------------------------------------------------

// Add a vertex on a polytope; the vertex is defined by 
// a domain coordinate (dc) and a range coordinate (rc).

vtkIdType vtkPolytopeGeometryForShapes::AddCoord(const double *dc, const double *rc)
{ 
  vtkIdType pid;
  int unseen;
  double *base;

  // The *domain* coordinate is unique to each polytope vertex;
  // check if the coordinate is already known, if so return the
  // existing coordinate id.
  unseen = this->ActiveDomCoord->InsertUniquePoint(dc, pid);
  // if this is a *new* domain point, store the corresponding range.
  if (unseen){
    for (int i = 0; i < N; i++)
      {
      this->ActiveScalars->InsertComponent(pid, i, rc[i]);
      }    
    }

  // Coordinate ids are represented as negative numbers.
  // Encode pid as a coordinate and return.
  return -1 - pid;
}

// ------------------------------------------------------------------


// Return a domain coordinate for an active vertex via a
// caller-allocated buffer.
void vtkPolytopeGeometryForShapes::GetActiveDomCoord(vtkIdType id, double *c)
{
  // pass buffer to PT object to fill.
  this->ActiveDomCoord->GetPoint(-1-id, c);
}

// ------------------------------------------------------------------

// Copy a range coordinate into caller-allocated buffer.
// TODO: consider removing need to copy by passing reference to 
//       coordinate data in VTK array.
void vtkPolytopeGeometryForShapes::GetActiveRngCoord(vtkIdType id, double *c)
{
  id = -1 - id;
  for (int i = 0; i < this->N; i++)
    {
    c[i] = this->ActiveScalars->GetComponent(id, i);
    }
}

// ------------------------------------------------------------------

// Return a specified component of an active vertex domain coordinate.
double vtkPolytopeGeometryForShapes::GetActiveDomCoordComponent(vtkIdType id,  int c)
{
  double coord[MaxDomainDim];
  this->GetActiveDomCoord(id, coord);
  return coord[c];
}

// ------------------------------------------------------------------

// Return a specified component of an active vertex range coordinate.
double vtkPolytopeGeometryForShapes::GetActiveRngCoordComponent(vtkIdType id,  int c)
{
  return this->ActiveScalars->GetComponent(id, c);
}

// ------------------------------------------------------------------

// Return a specified component of a STORED vertex range coordinate.
double vtkPolytopeGeometryForShapes::GetStoredRangeComponent(vtkIdType id, int c)
{
  return this->StoredScalars->GetComponent(id, c);
}

// ------------------------------------------------------------------

// Return the pre-computed lower corner of the bounding box for
// an active polytope.
double *vtkPolytopeGeometryForShapes::GetPtopeCoordMin(vtkIdType id)
{
  return this->Ptope[id].minCoord;
}

// ------------------------------------------------------------------

// Helper function:
// Compare two double values subject to an error tolerance.
// Return -1, 0, 1 respectively if t is less than,
// approximately equal to, or greater than, base.
int Compare_vtkPolytopeGeometryForShapes(double t, double base)
{
  if (fabs(t - base) < EPSILON)
    return 0;
  if (t < base)
    return -1;
  return 1;
}

// ------------------------------------------------------------------

// Access the facets of a given polytope by returning a pointer to
// the first facet and the number of facets in the polytope.
void vtkPolytopeGeometryForShapes::GetPtopeFacets(vtkIdType pid, vtkIdType *&ids, int &nr)
{
  ids = this->Facet + this->Ptope[pid].start;
  nr  = this->Ptope[pid].size;
}



// ------------------------------------------------------------------

// Return pointer to buffer holding the center of a given polytope;
// center is cached in the polytope, and is computed on demand. 
// Optimize special case of low-dimensional polytopes, but in general
// invoke function recursively on facets, averaging over facet centers.
double *vtkPolytopeGeometryForShapes::GetPtopeCenter(vtkIdType pid)
{
	//cout<<"\nthis->Ptope : "<<this->Ptope<<" , pid : "<<pid<<"\n";
  struct _Ptope *entry = this->Ptope + pid;
  double *fcenter, *base;
  double point[MaxDomainDim];
  vtkIdType a;
  vtkIdType *facet;

  if (entry->center)
    {
    // pre-computed point, return pointer to cache.
    return entry->center;
    }   
    
    // Changed to multi-dimensional point
    entry->center = new double[this->M];    
    
    
  // No cached center, so compute cache.
  // Initialize cache entry to zeroes.
  for (int i = 0; i < this->M; i++)
    {    
    entry->center[i] = 0.0;
    }

    
  // Branch on polytope dimension to handle common cases cheaply.    
  facet = this->Facet + this->Ptope[pid].start;
  switch (entry->dim) {

    case 0:
      // Point: center is just the point itself.
      a = -1 - facet[0];
      this->ActiveDomCoord->GetPoint(a, entry->center);
      break;
      
    case 1:
      // Edge: center is midway between endpoints.
      a = -1 - facet[0];
      this->ActiveDomCoord->GetPoint(a, point);
      a = -1 - facet[1];
      this->ActiveDomCoord->GetPoint(a, entry->center);
 
      for (int j = 0; j < this->M; j++)
        {
        entry->center[j] = 0.5 * (entry->center[j] + point[j]);
        }

      break;
      
    default:
      // Recursively compute facet centers and 
      // accumulate sum.
         
      for (int i = 0; i < entry->size; i++)
        {
        fcenter = this->GetPtopeCenter(facet[i]);
        for (int j = 0; j < this->M; j++)
          {
          entry->center[j] += fcenter[j];
          }
        }
      
      // Take average of facet centers
      // as polytope center.
      for (int j = 0; j < this->M; j++)
        {
        entry->center[j] /= entry->size;
        }
    }
  return entry->center;
}

// ------------------------------------------------------------------
// JCN-specific code:
// ------------------------------------------------------------------

// Split an active polytope against a sequence of increasing 
// thresholds, returning new polytope fragments and plane-polytope 
// intersections ("ends") via caller-allocated buffers.  Note
// that number of returned components will be one more than number
// of planes, due to a possible residue beyond the last plane.
// "nrts" refers to the number of planes, and does not include this
// extra entry.

void vtkPolytopeGeometryForShapes::Split(
  vtkIdType ptid,         // polytope id to split
  double *thresholds,     // ascending sequence of thresholds.
  int nrts,               // number of thresholds in buffer
  int field,              // component of domain being sliced
  vtkIdType *frags,       // return buffer for fragments
  vtkIdType *ends)        // return buffer for intersections
{
  struct _Ptope *entry;
  int t_comp, prev_comp;
  vtkIdType a;
  double ai;

  bool placed = false;    // has the ptope being fully processed?

  entry = this->Ptope + ptid;
  
  // Easy case: reuse cached result.
  if (entry->fragments && entry->ends && entry->splitDim == field)  
    {
    // This ptope has already been fragmented on this axis, 
    // so we re-use the fragment and end information cached with the ptope entry.
    for (int i = 0; i < nrts+1; i++)
      {
      frags[i] = entry->fragments[i];
      ends[i] = entry->ends[i];
      }
    return;
    }

  // Clear and reallocate cache: nr of thresholds will vary between polytopes
  // and axes.
  // NOTE: in principle could trade space for speed and have a fixed size
  // cache in every ptope entry - initial evidence was that speed gain is
  // minimal.

  if (entry->fragments)
    {
    delete [] entry->fragments;
    delete [] entry->ends;
    }

  // Allocate space for one fragment per cutting plane PLUS residue.

  entry->fragments = new vtkIdType[nrts + 1];
  entry->ends = new vtkIdType[nrts + 1];
  entry->splitDim = field;

  // Handle polytope based on dimension.

  switch (entry->dim) {

    case 0: 

      // Rare: polytope is a single point.
      // Point either:
      // - lies on a plane
      // - lies between two planes
      // - comes after the last plane.

      a = this->Facet[entry->start];
      ai = this->GetActiveRngCoordComponent(a, field);

      prev_comp = -1; // comparison result for previous plane. 
      for (int c = 0; c < nrts; c++)
        {
        // Compare coordinate component with current plane threshold.
        t_comp = Compare_vtkPolytopeGeometryForShapes(ai, thresholds[c]);
        switch (t_comp) {
          case -1: 
            // plane is > component.
            entry->fragments[c] = entry->ends[c] = ends[c] = frags[c] = 0; 
            break;
          case  0: 
            // component matches plane, point lies on plane.
            entry->ends[c] = ends[c] = a; 
            placed = true;
            entry->fragments[c] = frags[c] = 0; 
            break;
          case  1:
            // component is < plane 
            if (prev_comp < 0)
              {
              // previous plane was < component, so point
              // lies between both planes.
              entry->ends[c] = ends[c] = 0; 
              entry->fragments[c] = frags[c] = a;
              placed = true;
              }
            else
              {
              entry->ends[c] = entry->fragments[c] = ends[c] = frags[c] = 0;
              }
            break;
          }
        prev_comp = t_comp;
        }
      // handle case when point lies beyond last cutting plane.
      entry->ends[nrts] = ends[nrts] 
          = entry->fragments[nrts] = frags[nrts]
          = placed ? 0 : a;
      // END OF ZERO-D (point) CASE
      break;
      
    case 1: // EDGE 
      this->SplitEdge(entry, thresholds, nrts, field, frags, ends);
      break;

    default: // arbitrary dimension polytope.
      this->SplitPtope(entry, thresholds, nrts, field, frags, ends);

    } // switch dimension

} // Split.

// ------------------------------------------------------------------

// Split a polytope of dimension d > 1.
// Recursively split the polytope's facets, then
// at each threshold,
// - attempt to build a new "end" facet from the fragments
//   of the facet intersections
// - attempt to build a new fragment from:
//   1.  the fragments of the facets
//   2.  the end facet, if it exists
//   3.  the PREVIOUS end facet, if it exists.

void vtkPolytopeGeometryForShapes::SplitPtope(
  struct _Ptope *entry,   // the ptope to be split
  double *thresholds,     // cutting planes, ascending
  int nrts,               // number of cutting planes
  int field,              // component being cut
  vtkIdType *frags,       // buffer to hold fragments
  vtkIdType *ends)        // buffer to hold intersections.
{

  // Allocate temp space for storing output fragments and ends formed by 
  // cutting each part of this ptope against the thresholds.
  // There are nrts + 1 possible sets of output, as some parts of the
  // ptope may lie beyond the final threshold.
  // There are entrysize + 2 rows, as in addition to fragments lying 
  // between cutting planes, there are up to two end ptopes to consider.
  vtkIdType *rec_frags = new vtkIdType[(nrts + 1) * (entry->size + 2)];
  vtkIdType *rec_ends  = new vtkIdType[(nrts + 1) * (entry->size + 2)];

  // ids for the polytopes formed by fragments lying in cutting planes.
  vtkIdType end_tope, end_prev;
  int nrParts;
  
  // Recursively split each facet on the same field/thresholds.
  // Note offsets to start of buffer for each threshold include
  // the slots for the residues.
  for (int f = 0; f < entry->size; f++)
    {
    this->Split(
        this->Facet[entry->start + f],
        thresholds,      
        nrts,
        field,
        rec_frags + f*(nrts + 1),
        rec_ends  + f*(nrts + 1));
    } // for each face
  

  // Try to compute a new polytope formed by each cutting 
  // plane and by the residue (hence <= condition for loop 
  // termination).
  // We convert end fragments into a new polytope, adding
  // these to each pair of polytopes, using the next
  // two rows of the rec_ table ... so assumption is that
  // input polytopes are at least to facets smaller than
  // MaxNrFacetsPerPtope.

  end_prev = 0;
  for (int i = 0; i <= nrts; i++)
    {
    // Track number of *candidate* parts for each new polytope.
    // Each facet will contribute a candidate.
    nrParts = entry->size;    
    // Attempt to build an intersection polytope from end fragments
    // lying in the plane (for the residue this will always be nil).
    end_tope = entry->ends[i] = ends[i] 
        = this->Select(rec_ends, nrts+1, nrParts, i, entry->dim-1, entry->dim);
    // If there is an intersection polytope, add it to the
    // candidate components for making the next polytope fragment.
    if (end_tope)
      {
      rec_frags[(nrts + 1) * nrParts + i] = end_tope;
      nrParts++;
      }
    // If the *previous* cut resulted in an end-polytope, add it
    // to the candidate fragments for this polytope.
    if (end_prev)
      {
      rec_frags[(nrts + 1) * nrParts + i] = end_prev;  
      nrParts++;
      }
    // The end polytope from this plane will be the previous for the next plane.
    end_prev = end_tope;
    // Attempt to build a polytope fragment from the facet fragments and any 
    // end polytopes.
    entry->fragments[i] = frags[i] 
        = this->Select(rec_frags, nrts+1, nrParts, i, entry->dim, entry->dim);
    }

  // Recover temporary space.
  delete [] rec_frags;
  delete [] rec_ends;

} // SplitPtope

// ------------------------------------------------------------------

// Split an edge against a set of cutting planes to produce a 
// collection of segments.  Degenerate cases, e.g. where an edge
// lies within a plane or between planes, are handled.

void vtkPolytopeGeometryForShapes::SplitEdge(
  struct _Ptope *entry,   // The polytope to be split (will be an edge)
  double *thresholds,     // Thresholds on which to cut the edge
  int nrts,               // Number of thresholds provided
  int field,              // Coordinate component being cut.
  vtkIdType *frags,       // Buffer for fragments (includes a residue slot)
  vtkIdType *ends)        // Buffer for intersections
{
  double t;
  int t_comp, prev_comp;
  vtkIdType temp[2];      // temporary edge between two vertices.
  vtkIdType a, b;

  // Buffers for domain/range coordinates of edge endpoints.
  double domCoordA[MaxDomainDim];
  double domCoordB[MaxDomainDim];

  double rngCoordA[MaxRangeDim];
  double rngCoordB[MaxRangeDim];

  // Buffers for interpolated coordinates.
  double domCoord[MaxDomainDim];
  double rngCoord[MaxRangeDim];
  
  double ai, bi;
  vtkIdType coord_prev;
  bool working;

  // Get indices of edge endpoints, range coordinates of
  // endpoints, and values for component being cut against.
  a = this->Facet[entry->start];
  b = this->Facet[entry->start+1];
  this->GetActiveRngCoord(a, rngCoordA);
  this->GetActiveRngCoord(b, rngCoordB);
  ai = rngCoordA[field];
  bi = rngCoordB[field];

  // Ensure edge is oriented so that ai < bi
  if (bi < ai)
    {
    std::swap<double>(ai, bi);
    std::swap<vtkIdType>(a, b);
    for (int i = 0; i < this->N; i++)
      {
      std::swap<double>(rngCoordA[i], rngCoordB[i]);
      }
    }

  // Separate out the case of edge lying in the cutting plane.
  
  if (Compare_vtkPolytopeGeometryForShapes(ai, bi))
    { 
    // edge not coplanar, so seek intersections.  

    coord_prev = a;   // Last endpoint of edge to be processed
    working = true;   // Is output of edge ongoing?
    
    for (int c = 0; c < nrts; c++)
      {
      // Find parametric coordinate of edge/plane intersection,
      // and classify as below / on / above lowest endpoint.
      t = (thresholds[c] - ai) / (bi - ai);
      t_comp = Compare_vtkPolytopeGeometryForShapes(t, 0.0);

      switch (t_comp) {
        case -1: 
          // edge starts at component value above that of the plane.
          entry->fragments[c] = entry->ends[c] = frags[c] = ends[c] = 0;
          break;

        case  0: 
          // left-end of edge lies on the plane
          entry->fragments[c] = entry->ends[c] = frags[c] = ends[c] = 0; 
         break;                

        case  1: 
          // edge starts below plane.
          // consider intersection with right hand of edge.
          t_comp = Compare_vtkPolytopeGeometryForShapes(t, 1.0);             
          switch (t_comp) {
            case -1: 
              // Edge ends at component value above plane, so is cut.
              // Interpolate intersection points in both 
              // domain and range spaces.
              this->GetActiveDomCoord(a, domCoordA);
              this->GetActiveDomCoord(b, domCoordB);
           

              for (int j = 0; j < this->M+1; j++)
                {
                domCoord[j] = (1-t)*domCoordA[j] + t*domCoordB[j];
                }
              for (int j = 0; j < this->N; j++)
                {
                rngCoord[j] = (1-t)*rngCoordA[j] + t*rngCoordB[j];
                }           
              // Edge segment starts at the previous edge endpoint computed.         
              temp[0]  = coord_prev;
              // Edge segment ends at intersection point.  This intersection
              // point is also returned through the end/intersection buffer.
              entry->ends[c] = ends[c]  = temp[1] = this->AddCoord(domCoord, rngCoord);
              // New polytope is defined for the edge fragment.
              entry->fragments[c] = frags[c] = this->AddPtope(temp, 1, 2);
              // This intersection point now becomes the last edge endpoint seen.
              coord_prev = temp[1];
              break;

            case  0: 
              // edge ends at cutting plane
              // connect previous edge endpoint to right end of edge.
              temp[0]  = coord_prev;
              temp[1]  = b;
              entry->fragments[c] = frags[c] = this->AddPtope(temp, 1, 2);
              entry->ends[c] = ends[c]  = b;
              coord_prev = temp[1];
              // as edge endpoint has been seen, there will be no residue.
              working = false;
              break;
              
            case  1: 
              // edge ends at a component value below that of the 
              // current plane. If the end of the edge has not been seen,
              // there is still a segment to output.

              if (working) 
                {
                // connect previous endpoint to endpoint of edge.
                temp[0]  = coord_prev;
                temp[1]  = b;
                entry->fragments[c] = frags[c] = this->AddPtope(temp, 1, 2);
                working  = false;                    
                }
              else
                {
                entry->fragments[c] = frags[c] = 0;
                }
              entry->ends[c] = ends[c]  = 0;
              break;
            } // switch on right-end.

        } // switch on left-end.
        
      } // for each cutting plane

    // If we haven't reached the end of the line, there is a
    // residue beyond the last cutting plane.  
    if (working)
      {
      temp[0] = coord_prev;
      temp[1] = b;
      entry->fragments[nrts] = frags[nrts] = this->AddPtope(temp, 1, 2);
      // NB since line end didn't hit a plane, end is empty.
      entry->ends[nrts] = ends[nrts] = 0;
      }
    else
      {
      entry->fragments[nrts] = entry->ends[nrts] = frags[nrts] = ends[nrts] = 0;  
      }
      
    } // edge not parallel to planes
  else 
    { 
    // Edge lies in a plane parallel to the cutting planes for this component.
    // Scan through planes until the position of the edge is determined.
    prev_comp = -1;
    working = false;
    for (int c = 0; c < nrts; c++)
      {
      t_comp = Compare_vtkPolytopeGeometryForShapes(thresholds[c], ai);
      switch (t_comp) {
        case -1: 
          entry->ends[c] = entry->fragments[c] = ends[c] = frags[c] = 0; 
          break;
        case  0: 
          entry->ends[c] = ends[c] = (entry - this->Ptope);  // id of original ptope
          entry->fragments[c] = frags[c] = 0; 
          working = true;
          break;
        case  1: 
          if (prev_comp < 0)
            {
            entry->ends[c] = ends[c] = 0; 
            entry->fragments[c] = frags[c] = (entry - this->Ptope);
            working = true;
            }
          else
            {
            entry->ends[c] = entry->fragments[c] = ends[c] = frags[c] = 0;
            }
          break;
        }
      prev_comp = t_comp;
      }
    // by definition, residue edges parallel to cutting plane does not 
    // intersect any cutting plane.
    entry->ends[nrts] = ends[nrts] = 0;
    if (working)
      {
      entry->fragments[nrts] = frags[nrts] = 0;
      }
    else
      {
      entry->fragments[nrts] = frags[nrts] = (entry - this->Ptope); 
      }
    } // else edge parallel to cutting plane

} // SplitEdge.

// ------------------------------------------------------------------

// Identify a polytope of dimension between d and endDim, where
// endDim <= d, from a collection of polytope fragments held in 
// a 2D array pointed to by "pts", with "nrts" entries per row.
// Polytope fragments are to be taken from an array column 
// given by "slice".  Rationale for this data layout is to be
// be found in the implementation of "Split" given below.
// NOTE: the "nrts" parameter passed  to Select INCLUDES the
// extra slot allocated for the residue (see Split). 

vtkIdType vtkPolytopeGeometryForShapes::Select(
  vtkIdType *pts,
  int nrts, // NB here = number of planes + 1 for the residue
  int nrps,
  int slice,
  int d, int endDim)
{  
  struct _Ptope *entry;
  vtkIdType pt;
  int next, pos;
  vtkIdType temp[MaxNrSlicesPerPtope];
  int end = endDim > 0 ? endDim : 0;

  // Attempt to form the largest possible face.

  // For efficiency, separate cases for searching for points or
  // searching for higher-dim polytopes.
  // d = dimensionality of components sought.
  
  // EITHER we find a single ptope of the required dimension,
  // OR we find sufficient facets to build one.

  // 1. Attempt to find one ptope:
  next = 0;
  for (int i = 0; i < nrps; i++)
    {
    if ((pt = pts[i*nrts + slice]) >= 0 && this->Ptope[pt].dim == d)
      {
      temp[next++] = pt;
      }
    } // for each candidate component
  
  if (next == 1) { return temp[0]; }  // exactly one, existing, ptope found.
  
  // 2. Attempt to build from facets.
  // Separate case of 1D (lines from points) & 2D; as vertices are stored
  // as negative values their dimension is implicit, saving on inner-loop 
  // tests.

  next = 0;
  if (d > 1)
    { 
    // not a line, so must test candidate dimension explicitly.
    for (int i = 0; i < nrps; i++)
      {
      if ((pt = pts[i*nrts + slice]) > 0 && this->Ptope[pt].dim == d-1)
        {
        for (pos = 0; pos < next && temp[pos] != pt; pos++) {;} 
        if (pos == next)
          { 
          temp[next++] = pt;
          }
        }
      } // for each candidate point
    if (next > d) 
      { 
      // create a new ptope of one dimension up from its component faces.
      return this->AddPtope(temp, d, next); 
      }
    } // if not a line.
  else
    { 
    // working with a 1-D line, look for -ve components.
    for (int i = 0; i < nrps; i++)
      {
      if ((pt = pts[i*nrts + slice]) < 0)
        {
        for (pos = 0; pos < next && temp[pos] != pt; pos++) {;} 
        if (pos == next)
          { 
          temp[next++] = pt;
          }
        }
      } // for each candidate point
    if (next == 2)
      {
      return this->AddPtope(temp, 1, 2); 
      }
    } // line case.

  // The worst case: we have no suitable polytope;
  // return the degenerate Nil polytope.
  return 0;
}

// ------------------------------------------------------------------
// Copy polytope information from active to stored representation.
// "ids" contains the polytope id for the top-level of the polytopes
// to be copied, "nrPtopes" is the number of top-level ids to be copied.

void vtkPolytopeGeometryForShapes::StoreActivePolytopes(vtkIdType *ids, int nrPtopes)
{
//	cout<<"\nNumber of polytopes  : "<<nrPtopes<<"\n";
  vtkIdType pid;
  bool unseen;
  int nrLocalPoints = this->ActiveDomCoord->GetNumberOfPoints();
  struct _Ptope *entry, *subPtope;
  vtkIdType facet, facetId;
  std::queue<vtkIdType> ptopes = std::queue<vtkIdType>();  
  
  double *p;
  double div;
  double sc, minF, sw;
  
  vtkIdType *polyIds = NULL;
  int polyNr;
   
  // Process each top-level polytope in turn.
/*  for (int i = 0; i < nrPtopes; i++)
    {

    // entry refers to new top-level ptope.    
    entry = this->Ptope + ids[i];    
    p = this->GetPtopeCoordMin(ids[i]);
       
    // Copy polytope min coordinate into stored table.
    for (int j = 0; j < this->N; j++)
      {
        this->StoredScalars->InsertComponent(this->NrStoredPtopes + i, j, p[j]);
      }
     
    // For each facet of the polytope, obtain the facet center and insert into
    // table of stored centers.  Update CenterFacets table to link the polytope
    // to the facet center id; at most two polytopes share a facet, if this was
    // the first occurrence of the facet center point, initialize the second
    // polytope slot to -1, otherwise use the second slot
    //stops at i = 1 and j = 0
    for (int j = 0; j < entry->size; j++)
     {
      facet = this->Facet[entry->start + j];
      p = this->GetPtopeCenter(facet);
      unseen = this->CenterLocator->InsertUniquePoint(p, pid);
      subPtope =  this->Ptope+facet;      

      if (unseen)
        {
         // First occurrence of this center, use component 0.
         this->CenterFacets->InsertComponent(pid, 0, this->NrStoredPtopes+ i);
         this->CenterFacets->InsertComponent(pid, 1, -1);
         
         if (this->M==2) 
           {
            this->GetPtopeFacets(facet,polyIds,polyNr);  
            for (int m = 0; m <polyNr; m++)
               {
                this->CenterFacetsId->InsertComponent(pid, m, -1 - *(polyIds+m));
               }
           }
         else if (this->M==3)
           {
            this->CenterFacetsId->InsertComponent(pid, 0, facet);
           }
        }
      else
        {
         // Second occurrence, use component 1.
	double t = CenterFacets->GetComponent(pid, 1);
	if(t!=-1)
		cout<<"\ntest?????\n";
         this->CenterFacets->InsertComponent(pid, 1, this->NrStoredPtopes + i);
          
         if (this->M==2) 
           {
            this->GetPtopeFacets(facet,polyIds,polyNr);  
            for (int m = 0; m <polyNr; m++)
               {
                this->CenterFacetsId->InsertComponent(pid, m+2, -1 - *(polyIds+m));
               }
           }
         else if (this->M==3)
           {
            this->CenterFacetsId->InsertComponent(pid, 1, facet);
           }
        }
      }
      
    }*/
  for (int i = 0; i < nrPtopes; i++)
    {

	    // entry refers to new top-level ptope.    
	    entry = this->Ptope + ids[i];    
	    p = this->GetPtopeCoordMin(ids[i]);
	       
	    // Copy polytope min coordinate into stored table.
	    for (int j = 0; j < this->N; j++)
	      {
		this->StoredScalars->InsertComponent(this->NrStoredPtopes + i, j, p[j]);
	      }
	    
	    for (int j = 0; j < entry->size; j++)
	     {
	      facet = this->Facet[entry->start + j];
	      p = this->GetPtopeCenter(facet);
	      unseen = this->CenterLocator->InsertUniquePoint(p, pid);
	  //    subPtope =  this->Ptope+facet;      

	      if (unseen)
		{
			vector<int> triangles;
			triangles.push_back(this->NrStoredPtopes+ i);
			edgeTriangleMap[pid]=triangles;
		}
	      else
		{
			edgeTriangleMap[pid].push_back(this->NrStoredPtopes + i);
		}
	      }
      }
  this->NrStoredPtopes += nrPtopes;
} // StoreActivePolytopes.

// ------------------------------------------------------------------

// Implement standard VTK interface to print this object.
// Polytope structure is output in three forms:
// - vertex table, showing domain and range coordinates
// - polytope structure table, showing polytope facets
// - polytope vertex table, listing just the vertices of each polytope.

void vtkPolytopeGeometryForShapes::PrintSelf(ostream &os, vtkIndent indent)
{
  struct _Ptope *entry = this->Ptope;
  vtkIdType pt, f;
  double p[MaxDomainDim]; // CHANGE to allow N-d

  // A queue of polytope ids is used to unpack a polytope
  // down to the constituent points.

  std::queue<vtkIdType> ptopes = std::queue<vtkIdType>();  
  vtkSmartPointer<vtkIdList> coords = vtkSmartPointer<vtkIdList>::New();
  
  this->Superclass::PrintSelf(os,indent);
  os << "Dimensions: " << this->M << "->" << this->N << endl;
  os << "Ptope table -------------------------------" << endl;
  // Output polytope facet table.
  for (int i = 0; i < this->NextPtope; i++, entry++)
    {
    os << i << '\t' << entry->dim << '\t';
    f = entry->start;
    if (entry->dim > 1)
      {
      for (int j = 0; j < entry->size; j++)
        {
        os << this->Facet[f+j] << " ";
        }
      }
    else
      {
      for (int j = 0; j < entry->size; j++)
        {
        os << (-1 - this->Facet[f+j]) << " ";
        }
      }
    os << "\t\t";
    if (entry->minCoord)
      {
      for (int j = 0; j < this->N; j++)
        {
        os << entry->minCoord[j] << " ";
        }  
	    }
    else
      {
      os << "<no minCoord>";
      }
    if (entry->fragments)
      {
      os << "FRAGS" << endl;
      }
    else
      {
      os << "NO FRAGS" << endl;
      }
    }

  // Output polytope vertex table
  os << "Ptope vertices ----------------------------" << endl;
  for (int i = 0; i < this->NextPtope; i++)
    {
    
    // Perform a breadth-first traversal over polytope-facet
    // hierarchy.  
    ptopes.push(i);   // initialize queue to be top-level polytope
    coords->Reset();  // initialize set of points.

    while (!ptopes.empty())
      {
      // pop polytope to process from queue, and inspect
      // its dimensionality.
      pt = ptopes.front();
      ptopes.pop();
      entry = this->Ptope + pt;

      f = entry->start;
      if (entry->dim > 1)
        {
        // facet found, push onto queue for processing
        for (int j = 0; j < entry->size; j++)
          {
          ptopes.push(this->Facet[f+j]);
          }
        }
      else
        {
        // edge/point found, facets are coordinates which are
        // added to set of coordinates.
        for (int j = 0; j < entry->size; j++)
          {
          coords->InsertUniqueId(-1 - this->Facet[f+j]);
          }
        }
      }
    os << i << '\t' << this->Ptope[i].dim << '\t';
    for (int i = 0; i < coords->GetNumberOfIds(); i++)
      {
      os << coords->GetId(i) << " ";
      }    
    os << endl;  
    }
  
  // Output coordinate table.
  os << "Coord table -------------------------------" << endl;
  
  for (int i = 0; i < ActiveDomCoord->GetNumberOfPoints(); i++)
    {
    os << i << "\t";
    ActiveDomCoord->GetPoint(i, p);
    for (int j = 0; j < this->M; j++)
       {
       os << p[j] << " ";
       }
     os << "\t|\t";
     for (int j = 0; j < this->N; j++)
        {
        os << this->ActiveScalars->GetComponent(i, j) << " ";
        }
      os << endl;
    }
    os << "-------------------------------------------" << endl;
}

// ------------------------------------------------------------------

// Specialised version of printing, showing polytope vertices and
// vertex coordinate information for a specified set of polytopes
// passed as a parameter.

void vtkPolytopeGeometryForShapes::PrintToplevel(ostream &os, vtkIdType *top, int sz)
{
  struct _Ptope *entry = this->Ptope;
  vtkIdType pt, f;
  double p[MaxDomainDim]; // Change

  std::queue<vtkIdType> ptopes = std::queue<vtkIdType>();  
  vtkSmartPointer<vtkIdList> coords = vtkSmartPointer<vtkIdList>::New();
  
  
  os << "Ptope vertices ----------------------------" << endl;
  
  // For each specified polytope, perform a breadth-first
  // traversal through the facet hierarchy to obtain a set 
  // of vertex ids.
  for (int i = 0; i < sz; i++)
    {
    // Initialize traversal queue and vertex set.
    ptopes.push(top[i]);
    coords->Reset();

    while (!ptopes.empty())
      {
      // dequeue a polytope and inspect dimensionality.
      pt = ptopes.front();
      ptopes.pop();
      entry = this->Ptope + pt;
      f = entry->start;
      if (entry->dim > 1)
        {
        // Polytope: add facets to queue.
        for (int j = 0; j < entry->size; j++)
          {
          ptopes.push(this->Facet[f+j]);
          }
        }
      else
        {
        // Edge/point: add vertices to vertex set.
        for (int j = 0; j < entry->size; j++)
          {
          coords->InsertUniqueId(-1 - this->Facet[f+j]);
          }
        }
      }
    os << i <<"(" << top[i] << ")\t" << this->Ptope[top[i]].dim << '\t';
    for (int i = 0; i < coords->GetNumberOfIds(); i++)
      {
      os << coords->GetId(i) << " ";
      }    
    os << "\t\t";
    if (entry->minCoord)
      {
      for (int i = 0; i < this->N; i++)
        {
        os << entry->minCoord[i] << " ";
        }
      }  
    else
      {
      os << "<no minCoord>" << endl;
      }

     os << endl;  
    }  
  os << "Coord table -------------------------------" << endl;
 
  for (int i = 0; i < ActiveDomCoord->GetNumberOfPoints(); i++)
     {
     os << i << "\t";
     ActiveDomCoord->GetPoint(i, p);
     for (int j = 0; j < this->M; j++)
        {
        os << p[j] << " ";
        }
      os << "\t|\t";
      for (int j = 0; j < this->N; j++)
         {
         os << this->ActiveScalars->GetComponent(i, j) << " ";
         }
      os << endl;
     }
  os << "-------------------------------------------" << endl;
}



