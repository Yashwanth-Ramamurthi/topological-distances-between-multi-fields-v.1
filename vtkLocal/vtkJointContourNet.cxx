/*=========================================================================

 *	File: vtkJointContourNet.cxx
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
#include "vtkBitArray.h"
#include "vtkCellData.h"
#include "vtkCell.h"
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkExecutive.h"
#include "vtkIndent.h"
#include "vtkIntArray.h"
#include "vtkIdList.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIdTypeArray.h"
#include "vtkJointContourNet.h"
#include "vtkMutableUndirectedGraph.h"
#include "vtkMergePoints.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolytopeGeometry.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkShortArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUndirectedGraph.h"
#include <vector>
#include <set>
#include <ostream>

// Tolerance used when computing slab coordinates.
#define EPSILON 1.0e-10

vtkStandardNewMacro(vtkJointContourNet);

//----------------------------------------------------------------------------
// Field structure and interfaces.
//----------------------------------------------------------------------------

// type for recording field information, both attributes provided by the
// user through external API, and variables needed at runtime.

typedef struct _Field {
  _Field() 
  {
    name = NULL; 
    pointScalars = NULL;
    slabNumbers  = NULL;
    nodeScalars  = NULL;
    fieldBase = 0.0;
    useBase = false;
  }
  
  ~_Field() 
  { 
    if (name) 
      {
      delete [] name;
      name = NULL;
      }
    pointScalars = NULL;
    if (slabNumbers)
      {
      slabNumbers->Delete();
      }
    if (nodeScalars)
      {
      nodeScalars->Delete();
      }
  }
    
  char *name;                                     // field name, specified by user.
  vtkDataArray                   *pointScalars;   // scalar value for this field at each input vertex.
  vtkSmartPointer<vtkShortArray> slabNumbers;     // quantization level for slabs
  vtkSmartPointer<vtkDoubleArray> nodeScalars;    // quantized scalar value placed into output graph.
  double slabWidth;                               // quantization level for this field.
  double fieldBase;                               // base from which to compute slabs
  bool useBase;                                   // use specified base rather than array minimum.
  double range[2];                                // global value range (min, max) for this field.
  } 
  Field;

// ------------------------------------------------------------------

// Add a new field to be used for computing the JCN.
// Caller's responsibility to ensure field is not duplicated.
void vtkJointContourNet::AddField(char *nm, double wd)
{
  Field *f = new Field();
  f->name = new char[strlen(nm) + 1];
  strcpy(f->name, nm);
  f->slabWidth = wd;
  f->useBase = false;  
  this->Fields.push_back(f);
}

// ------------------------------------------------------------------

// Add a new field to be used for computing the JCN, setting the
// base for slab computation.
// Caller's responsibility to ensure field is not duplicated.
void vtkJointContourNet::AddField(char *nm, double wd, double base)
{
  Field *f = new Field();
  f->name = new char[strlen(nm) + 1];
  strcpy(f->name, nm);
  f->slabWidth = wd;
  f->useBase = true;  
  f->fieldBase = base;
  this->Fields.push_back(f);
}

// ------------------------------------------------------------------
// Return name of nth field.
char * vtkJointContourNet::GetField(int n)
{
  if (n < this->Fields.size())
    {
    return this->Fields[n]->name;
    }
  else
    {
    vtkWarningMacro("Attempting to access non-existent field");
    return NULL;
    }
}


// ------------------------------------------------------------------
// Return number of fields.
int vtkJointContourNet::GetNumberOfFields()
{
  return this->Fields.size();
}

// ------------------------------------------------------------------

// Clear all fields.
void vtkJointContourNet::ClearFields()
{
  this->Fields.clear();
}


// ------------------------------------------------------------------
// Filter construction/destruction
// ------------------------------------------------------------------

vtkJointContourNet::vtkJointContourNet()
{
//  this->GetInformation()->Set(vtkAlgorithm::PRESERVES_RANGES(), 1);
  this->CenterFirstSlab = false;
  this->Fields = std::vector<Field*>();

  this->UseCellTopologyArray = false;
  this->UseSpatialCoordinateArray = false;

  this->TopologyArrayName = NULL;
  this->SpatialArrayName = NULL;

  this->IdentifyBoundarySlabs = false;
  this->M = -1;
  this->N = -1;
  this->UF = vtkSmartPointer<vtkIdTypeArray>::New();
  this->SetNumberOfOutputPorts(3);
}

// ------------------------------------------------------------------

vtkJointContourNet::~vtkJointContourNet()
{
  while (this->Fields.size() > 0)
    {
    this->Fields.pop_back();
    }
  this->Fields.clear();
  this->SetTopologyArrayName(0);
  this->SetSpatialArrayName(0);
}

// ------------------------------------------------------------------

// Implement VTK interface for displaying filter parameters.

void vtkJointContourNet::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << "Spatial (domain) dimension: " << this->M << endl;
  os << "Data (range) dimension: " << this->N << endl;
  os << indent << "Field names: ";
  if (this->Fields.size() < 1) 
    {
    os << "none defined." << endl;
    }
  else
    {
    os << endl;
    for (std::vector<Field*>::iterator it = this->Fields.begin();
        it != this->Fields.end();
        it++)
      {
      os << indent.GetNextIndent() << ((*it)->name);
      os << ", width " << (*it)->slabWidth << endl;
      }
    }
  os << "Use cell topology array: " << this->UseCellTopologyArray << endl;
  os << "Use spatial coord array: " << this->UseSpatialCoordinateArray << endl;
  os << "Cell topology array: " << (this->TopologyArrayName ? this->TopologyArrayName : "(null)") << endl;
  os << "Spatial coord array: " << (this->SpatialArrayName ? this->SpatialArrayName : "(null)") << endl;
  os << "Center first slab: " << this->CenterFirstSlab << endl;
}

// ------------------------------------------------------------------
// Union-find
// ------------------------------------------------------------------

// We use Tarjan's union-find algorithm.  The partition is implemented
// as an array.

static vtkSmartPointer<vtkIdTypeArray> UF;

// Find the root of a value stored in the union-find
// structure, performing path compression as we go.
vtkIdType vtkJointContourNet::Find(vtkIdType x)
{
  vtkIdType entry = this->UF->GetValue(x);
  if (entry < 0)
    {
    return x;
    }
  else
    {
    vtkIdType root = this->Find(entry);
    this->UF->SetValue(x, root);
    return root;
    }
}


//----------------------------------------------------------------------------
// VTK Pipeline Management
// Methods in the following section implement pipeline requests, and
// follow standard VTK patterns.
//----------------------------------------------------------------------------

int vtkJointContourNet::FillInputPortInformation(int, vtkInformation* info)
{
  // This filter uses the vtkDataSet cell traversal methods so it
  // suppors any data set type as input.
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
  return 1;
}

// ------------------------------------------------------------------

int vtkJointContourNet::FillOutputPortInformation(int ignored, vtkInformation* info)
{
  // Port 0 = JCN graph
  // Port 1 = optional unstructured grid of fragments.
  // Port 2 = optional unstructured grid of slabs boundary.
  if (ignored == 0)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUndirectedGraph");
    return 1;
    }
  if (ignored == 1)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
    }
   if (ignored == 2)
    {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
    }
  return 1;
}

// ------------------------------------------------------------------

vtkDataObject* vtkJointContourNet::GetOutput(int index)
{
  vtkDataObject *obj;
  if (index == 0)
    {
    obj = vtkUndirectedGraph::SafeDownCast(this->GetOutputDataObject(0));
    return obj;
    }
  if (index == 1)
    {
    obj = vtkUnstructuredGrid::SafeDownCast(this->GetOutputDataObject(1));
    return obj;
    } 
   if (index == 2)
    {
    obj = vtkUnstructuredGrid::SafeDownCast(this->GetOutputDataObject(2));
    return obj;
    } 
  return NULL;
}

// ------------------------------------------------------------------

vtkUndirectedGraph *vtkJointContourNet::GetTopologyGraph()
{
  return vtkUndirectedGraph::SafeDownCast(this->GetOutputDataObject(0));
}

// ------------------------------------------------------------------

int vtkJointContourNet::RequestDataObject(
  vtkInformation* request,
  vtkInformationVector** inputVector ,
  vtkInformationVector* outputVector)
{
  vtkInformation* info0 = outputVector->GetInformationObject(0);
  vtkInformation* info1 = outputVector->GetInformationObject(1);
  vtkInformation* info2 = outputVector->GetInformationObject(2);
  vtkUndirectedGraph *output = vtkUndirectedGraph::SafeDownCast(
        info0->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *outputFrags = vtkUnstructuredGrid::SafeDownCast(
        info1->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *outputBoundary = vtkUnstructuredGrid::SafeDownCast(
        info2->Get(vtkDataObject::DATA_OBJECT()));

  if (!output)
    {
    output = vtkUndirectedGraph::New();
    info0->Set(vtkDataObject::DATA_OBJECT(), output);
    output->Delete();
    }

  if (!outputFrags)
    {
    outputFrags = vtkUnstructuredGrid::New();
    info1->Set(vtkDataObject::DATA_OBJECT(),outputFrags);
    outputFrags->Delete();
    }
    
  if (!outputBoundary)
    {
    outputBoundary = vtkUnstructuredGrid::New();
    info2->Set(vtkDataObject::DATA_OBJECT(), outputBoundary);
    outputBoundary->Delete();
    }

  return 1;
}



// ------------------------------------------------------------------
// The JCN Algorithm.
// ------------------------------------------------------------------

// Convert a value v in field f into the corresponding slab coordinate
// by applying the level of quantization required for that field.
// If used, first slab centering is assumed to have already been 
// accounted for in the field minimum.

double vtkJointContourNet::SlabCoordinate(int f, double v)
{
  double minF = this->Fields[f]->range[0];
  double sw   = this->Fields[f]->slabWidth;
  double div  = (v + EPSILON - minF) / sw;

  return minF + sw*floor(div);
}


// ------------------------------------------------------------------
// Compute the Joint Contour Net.
// This method makes extensive use of vtkPolytopeGeometry, and
// familiarity with that call is assumed.
//
// Overview:
// ------------------------------------------------------------------
// For each input cell C
//   construct a poplytope [simplex] P from C
//   initialise a working set S of polytope fragments to {P}
//   for each field f
//      fragment each polytope in S
//   place fragments into polytope storage
//
// Initialize union-find structure with each fragment in 
//     a separate class.
// For each fragment boundary:
//   if the neighbouring fragments are equivalent
//     merge fragment partitions in the union-find
// Count number of partitions (== number of slabs)
// Construct graph
//   one node per slab
//   one edge between neighbouring fragments in different partitions.
// ------------------------------------------------------------------

int vtkJointContourNet::RequestData(vtkInformation*,
                                 vtkInformationVector** inputVector,
                                 vtkInformationVector* outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkInformation *outInfo1 = outputVector->GetInformationObject(1);
  vtkInformation *outInfo2 = outputVector->GetInformationObject(2);
  
  // ------------------------------------------------------
  // PIPELINE: input and output data sets
  //           point & cell data, and points.
  //           multifield dimensionality
  // ------------------------------------------------------
  
  vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(
      inInfo->Get(vtkDataObject::DATA_OBJECT()));
      
  vtkGraph *output = vtkGraph::SafeDownCast(
      outInfo->Get(vtkDataObject::DATA_OBJECT()));
      
  vtkUnstructuredGrid *outFrags = vtkUnstructuredGrid::SafeDownCast(
        outInfo1->Get(vtkDataObject::DATA_OBJECT()));  
        
  vtkUnstructuredGrid *outBoundary = vtkUnstructuredGrid::SafeDownCast(
        outInfo2->Get(vtkDataObject::DATA_OBJECT()));    

  vtkDataSetAttributes *inPD     = input->GetPointData();
  
  vtkSmartPointer<vtkUnstructuredGrid> cIdList = 
                  vtkSmartPointer<vtkUnstructuredGrid>::New();

  
  vtkPointData *outFragPD, *outBndrPD;
  vtkCellData  *outFragCD, *outBndrCD;
    
  // If ExplicitTopology turned on, then we read the input from
  // a array named "Cells". This array can be set using
  // vtkSetFieldData method. 
  vtkDataSetAttributes *inCD;

  if (this->UseCellTopologyArray)
    {
    inCD = vtkDataSetAttributes::SafeDownCast(input->GetFieldData());
    } 
  else 
    {
    inCD = input->GetCellData();
    }
    
  // Attempt to find a cell-data array called Cells, in case
  // explicit topology is to be used.
  vtkIdTypeArray *cells = NULL;

  vtkIdType numPts;
  vtkIdType numCells;

  if (this->UseCellTopologyArray)
    {
    if (!this->TopologyArrayName)
      {
      vtkErrorMacro("No array name specified for cell topology.");
      return 0;
      }
    cells = vtkIdTypeArray::SafeDownCast(
        inCD->GetArray(this->TopologyArrayName));
    if (!cells)
      {
      vtkErrorMacro("No topology cell array on input.");
      return 0;
      }
    numCells = cells->GetNumberOfTuples();
    this->M  = cells->GetNumberOfComponents() - 1;
    }
  else
    {
    numCells = input->GetNumberOfCells();
    this->M  = input->GetCell(0)->GetNumberOfPoints() - 1;
    }
  this->N = this->Fields.size();
  
  if (this->M > 3 && !this->UseSpatialCoordinateArray)
    {
    vtkErrorMacro("Spatial dimension > 3 requires separate coordinate array.");
    return 0;
    }

  vtkDataArray *points = NULL;
  if (this->UseSpatialCoordinateArray)
    {
    if (!this->SpatialArrayName)
      {
      vtkErrorMacro("Spatial array name not specified.");
      return 0;
      }
    points = input->GetPointData()->GetArray(this->SpatialArrayName);
    if (!points)
      {
      vtkErrorMacro("Spatial array not found.");
      return 0;
      }
    }
  else
    {
    points = input->GetPoints()->GetData();
    }
  numPts = points->GetNumberOfTuples();
  
  vtkIdType estOutputSize = numCells * 5;
  
  outFrags->Allocate(estOutputSize);
  cIdList->Allocate(estOutputSize);
  outFragPD = outFrags->GetPointData();
  outFragCD = outFrags->GetCellData();
  
  outBoundary->Allocate(estOutputSize);
  outBndrPD = outBoundary->GetPointData();
  outBndrCD = outBoundary->GetCellData();
  
  double dataBounds[MaxDomainDim];

  // Find the range of each spatial dimension
  // Ranges are stored as consecutive min/max pairs in
  // the bounds array.
  for (int d = 0; d < this->M; d++)
     {
      points->GetRange(dataBounds + d*2, d);
     }

   
  // ------------------------------------------------------
  // Check for pathological conditions, and initialize
  // output structures.
  // ------------------------------------------------------

  if (numCells < 1 || numPts < 1)
    {
    vtkErrorMacro("No sample points in input data.");
    return 0;
    }
  if (this->N < 1) 
    {
    vtkErrorMacro("No data fields.");
    return 0;
    }

  // ------------------------------------------------------
  // Polytope storage.
  // ------------------------------------------------------

  vtkSmartPointer<vtkPolytopeGeometry> geometry 
      = vtkSmartPointer<vtkPolytopeGeometry>::New();

  // ------------------------------------------------------
  // Graph construction: data structures and support.
  // ------------------------------------------------------

  vtkSmartPointer<vtkMutableUndirectedGraph> topograph
      = vtkSmartPointer<vtkMutableUndirectedGraph>::New();

  vtkSmartPointer<vtkIdTypeArray> cellVertex
      = vtkSmartPointer<vtkIdTypeArray>::New();
      
   // Arrays for tracking number of fragments in each slab,
  // then providing this information as vertex data.

  // Tracks the number of facets in a fragment. 
  vtkSmartPointer<vtkIntArray> fragSize
      = vtkSmartPointer<vtkIntArray>::New();

  // Tracks the number of fragments in a slab. 
  vtkSmartPointer<vtkIntArray> nodeSize
      = vtkSmartPointer<vtkIntArray>::New();
      
  // Tracks the number of facets in a slab.     
  vtkSmartPointer<vtkIntArray> slabSize
      = vtkSmartPointer<vtkIntArray>::New();
      
  // Center of the frags    
  vtkSmartPointer<vtkDoubleArray> fragCenters 
      = vtkSmartPointer<vtkDoubleArray>::New();


 // Points for the output grid.  Output points will be a superset of the
  // input, so allocate and deep copy the input points.  Output will also
  // have the same bounds as the input.
  vtkIdType ignored;
  vtkPoints *oldPoints = input->GetPoints();
  vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
  
  // Use floats for points to reduce memory o/head:
  // may need to review in case of artefacts.
  newPoints->Allocate(numCells*5, numPts/2);

  vtkSmartPointer<vtkMergePoints> pointLocator = vtkSmartPointer<vtkMergePoints>::New();
  pointLocator->InitPointInsertion (newPoints, input->GetBounds());


  
  for (int i = 0; i < this->N; i++)
    {
       vtkSmartPointer<vtkDoubleArray> arry = vtkSmartPointer<vtkDoubleArray>::New();
       arry->SetName(this->Fields[i]->name);
       outFragPD->AddArray(arry);
       outBndrPD->AddArray(arry);
    }

 
  // ------------------------------------------------------
  // Polytope information.
  // Array dimensions use constants from vtkPolytopeGeometry.
  // ------------------------------------------------------

  double domRange[MaxDomainDim*2];
  double rngRange[MaxRangeDim*2];

  vtkIdType cellPoints[MaxDomainDim+1];
  double coord[MaxDomainDim][MaxRangeDim];
  double fv;
  double sw;
  double minF;
  double thresholds[MaxNrSlicesPerPtope];
  double point[MaxDomainDim];  //Dimension M
  
  int slice;
  
  // Buffers for holding fragment and intersection
  // polytopes returned by vtkPolytopeGeometry::Slice.

  vtkIdType frags[MaxNrSlicesPerPtope];
  vtkIdType ends[MaxNrSlicesPerPtope];

  // Buffers for holding sets of top-level polytopes.  
  // One buffer holds the input polytopes for fragmenting
  // against the next field, the other buffer collects the
  // output fragments.
  // Pointers are used to swapover buffer roles between fields.
  vtkIdType top1[MaxNrSlicesPerPtope];
  vtkIdType top2[MaxNrSlicesPerPtope];
  vtkIdType *toFragment   = top1;
  vtkIdType *newFragments = top2;

  int nrToFragment   = 0;
  int nrNewFragments = 0;
  
  vtkIdType pointIds[MaxDomainDim];
  vtkSmartPointer<vtkIdList> pids = vtkSmartPointer<vtkIdList>::New();
  vtkIdType pID, base, end;
  double *pcoord; // polytope minimum coordinate.

  // ------------------------------------------------------
  // Traverse named scalar fields and initialize remainder of
  // structure from input dataset.
  // ------------------------------------------------------
  
  vtkSmartPointer<vtkDataArray> array;
  for (std::vector<Field *>::iterator it = this->Fields.begin();
      it != this->Fields.end();
      it++)
    {
    array = inPD->GetArray((*it)->name);
    if (array)
      {
      (*it)->pointScalars = array;
      array->Modified();
      array->GetRange((*it)->range);
      // Allocate storage for slab number and node scalar arrays.
      (*it)->slabNumbers  = vtkSmartPointer<vtkShortArray>::New();
      (*it)->nodeScalars  = vtkSmartPointer<vtkDoubleArray>::New();
      (*it)->nodeScalars->SetName(array->GetName());

      // Handle first-slab centering: if the first slab is centered 
      // over the lowest sample in the field, the effect is the same
      // as slicing a field whose minimum value is a half slab-width
      // lower.
      if ((*it)->useBase)
        {
        (*it)->range[0] = (*it)->fieldBase;
        }
      if (this->CenterFirstSlab)
        {
        (*it)->range[0] -= (*it)->slabWidth/2;
        }
      }
    else
      {
      vtkErrorMacro("Array not found in input point data: " << (*it)->name);
      return 0;
      }
    }
    
  // Initialize the geometry processor with dimensionality and 
  // global spatial bounds.

  geometry->Initialize(this->M, this->N, dataBounds);

  // ------------------------------------------------------
  // MAIN LOOP:
  //
  // FOR EACH cell c in the dataset ...
  // ------------------------------------------------------
  for (vtkIdType c = 0; c < numCells; c++)
    {
    // Obtain domain coordinate for first point in the cell.
    // Source of point id depends on choice between explicit
    // or implicit topology.
    if (this->UseCellTopologyArray)
      {
      points->GetTuple(cells->GetComponent(c, 0), point);
      }
    else
      {
      input->GetCellPoints(c, pids);
      points->GetTuple(pids->GetId(0), point);
      }
    // Compute the domain bounds of the current simplex.
    // Initialize range min/max to first point ordinates.

    for (int i = 0; i < this->M; i++)
      {
      domRange[i*2] = domRange[i*2+1] = point[i];    
      }
    // Compare bounds with remaining M points, adjusting
    // bounds as necessary.      
    for (int j = 1; j < this->M+1; j++)
      {
      if (this->UseCellTopologyArray)
        {
        points->GetTuple(cells->GetComponent(c, j), point);
        }
      else
        {
        points->GetTuple(pids->GetId(j), point);
        }
      for (int i = 0; i < this->M; i++)
        {
        domRange[i*2]   = std::min(domRange[i*2  ], point[i]);
        domRange[i*2+1] = std::max(domRange[i*2+1], point[i]);
        }
      }
    // Prepare the geometry processor to receive new active polytopes.
    geometry->ResetForNextCell(domRange);
        
    // Collect the samples assigned to each point in the simplex,
    // computing the range in each range dimension, and creating 
    // the corresponding point in the geometry processor.      
      
    for (int i = 0; i < this->N; i++)
      {
      rngRange[i*2]   = VTK_DOUBLE_MAX;
      rngRange[i*2+1] = VTK_DOUBLE_MIN;
      }
    for (int v = 0; v < this->M + 1; v++)
      {  
      // Initialize the simplex with its point ids. In the subsequent 
      // step, the (M+1) (dim-1)-dimensional facets will be constructed 
      // as simplicial subsets of these points.
      pID = this->UseCellTopologyArray
                ? cells->GetComponent(c, v) 
                : pids->GetId(v);
 
      points->GetTuple(pID, point);
     
      // Initialize the simplex with its point ids. In the subsequent 
      // step, the (M+1) (dim-1)-dimensional facets will be constructed 
      // as simplicial subsets of these points.

      // Cache values for the cell, and compute min/max for each
      // field within the cell.
      // Initial point position
    
      for (int i = 0; i < this->N; i++)
        {
        coord[v][i]
            = fv = this->Fields[i]->pointScalars->GetComponent(pID, 0);
        if (fv < rngRange[i*2  ]) rngRange[i*2  ] = fv;
        if (fv > rngRange[i*2+1]) rngRange[i*2+1] = fv;
        }
       // Point index for each cell // Add domain point and range points
       pointIds[v] = geometry->AddCoord(point, coord[v]);
       
      }  // for each domain point v in each cell, add domain and range values by AddCoord
      
    // Create a simplex from the cell points; this single polytope
    // forms the initial input to fragmentation.
    nrToFragment = 0;
    toFragment[nrToFragment++] 
        = geometry->BuildFromSimplex(pointIds, this->M);
         
    // For each field, compute the threshold values at which the
    // field intersects the cell, then slice the input polytopes 
    // against those thresholds.    
    for (int f = 0; f < this->Fields.size(); f++)
      {
      
        sw = this->Fields[f]->slabWidth;
      
        // First threshold is largest slab that contains mimumum
        // value of this component for the cell.
        thresholds[0] = this->SlabCoordinate(f, rngRange[f*2]);     
         
        // Attempt to add further thresholds, until last threshold
        // lies beyond largest component value in cell.
        slice = 1;      
        for (; thresholds[slice-1] < rngRange[f*2+1]; slice++)
           {
              thresholds[slice] = thresholds[slice-1] + sw;
           }

      // Trivial test: the polytope ISN'T cut, as the only cutting plane
      // lies at or before the cell minimum (i.e. the cell is fully
      // contained within one quantization interval).
      if (slice == 1 && thresholds[0] <= rngRange[f*2])
        {
        // Trivial case: copy all input polytopes to output buffer.
        nrNewFragments = nrToFragment;
        for (int p = 0; p < nrToFragment; p++)
          {
          newFragments[p] = toFragment[p];
          }
        }
      else
        {
        // Non-trivial case:
        // Fragment each input polytope, adding output fragments
        // to output buffer.
        nrNewFragments = 0;
        for (int p = 0; p < nrToFragment; p++)
          {
          geometry->Split( 
              toFragment[p],      // id of ptope to be split
              thresholds,         // thresholds on which to split
              slice,              // number of thresholds
              f,                  // id of field used to split
              frags,              // buffer to hold fragment ids
              ends);              // buffer for ids of end-faces

          // Copy over fragments, including any residue (hence <= slice)
          for (int t = 0; t <= slice; t++)
            {
            if (frags[t])
              {
                 newFragments[nrNewFragments++] = frags[t];
              }
            }
          } // for each input fragment.
           
        } // non-trivial case.
         
      // newFragments is now populated with list of fragments.
      // swap buffers so that newFrags become the input for
      // fragmenting on the next field.
      std::swap<vtkIdType*>(toFragment, newFragments);
      std::swap<int>(nrToFragment, nrNewFragments);
    
      } // for each field
    // Copy all fragments for this cell from active to stored
    // polytopes.
    geometry->StoreActivePolytopes(toFragment, nrToFragment);
    

    vtkIdType *polyIds = NULL, 
              *faceIds = NULL, 
              *lineIds = NULL;   
    // Retrieve and store the fragment geometry.                    
    for (vtkIdType co = 0; co < nrToFragment; co++)
      {    
        vtkIdType pid = toFragment[co];
        int polyNr = 0;
        geometry->GetPtopeFacets(pid,polyIds,polyNr);  

        if (this->M == 2)
          {
           vtkSmartPointer<vtkIdList> cId = vtkSmartPointer<vtkIdList>::New();  
           // 2D case
           vtkSmartPointer<vtkIdList> face = vtkSmartPointer<vtkIdList>::New();         
           vtkIdType list[polyNr];
           vtkIdType replace[polyNr*2];  
           // number of vertex ids in a face
           for (int id=0; id < polyNr; id++)
              {            
                int faceNr = 0;
                geometry->GetPtopeFacets(*(polyIds+id),faceIds,faceNr);
                int temp1,temp2;
                //for each line segment
                for (int pnr = 0; pnr < faceNr; pnr++)
                   {   
                   replace[id*2+pnr] = -1 - *(faceIds+pnr);       
                   }                                
              }             
           list[0] = replace[0];
           list[1] = replace[1];
           replace[0] = -1;
           replace[1] = -1; 
           // Sort the vertices into right order for a polytope
           for (int q = 2; q < polyNr; q++)
              {
              int value = list[q-1];
              for (int p = 1; p < polyNr; p++)
                 {
                 int tempStart = replace[p*2];
                 int tempEnd = replace[p*2+1];
                 if (tempStart > 0)
                   {
                   if (tempStart == value) 
                     { 
                      list[q] = tempEnd;
                      replace[p*2]=replace[p*2+1] = -1; 
                     }
                   else if (tempEnd == value) 
                     {
                      list[q] = tempStart;
                      replace[p*2]=replace[p*2+1] = -1;        
                      }
                    }
                  }                     
              } //end vertex order
                 
             for (int q = 0; q < polyNr; q++)
                {
                  vtkIdType pid;
                  double p [] = {0,0,0};
                  double rangeP[this->N];
                  geometry->GetActiveDomCoord(-1-list[q], p);
                  geometry->GetActiveRngCoord(-1-list[q], rangeP);
                  // if it is a new point
                  if (pointLocator->InsertUniquePoint(p,pid))
                    {
                     for (int i = 0; i < outFragPD->GetNumberOfArrays(); i++)
                        {
                         outFragPD->GetArray(i)->InsertNextTuple1(rangeP[i]);   
                        }
                    }
                  int num;
                  face->InsertNextId(pid);
                  cId->InsertNextId(list[q]);   
                }            
            outFrags->InsertNextCell(VTK_POLYGON, face); 
            cIdList->InsertNextCell(VTK_POLYGON,cId);
        
          } // 2D case
        else if (this->M == 3) 
          { 
            // 3D case
            vtkSmartPointer<vtkIdList> poly = vtkSmartPointer<vtkIdList>::New();
            vtkSmartPointer<vtkIdList> polyId = vtkSmartPointer<vtkIdList>::New();
            // number of faces
            poly->InsertNextId(polyNr); 
   
           for (int id=0; id < polyNr; id++)
              {            
                int faceNr = 0;
                geometry->GetPtopeFacets(*(polyIds+id),faceIds,faceNr);
                poly->InsertNextId(faceNr);
                polyId->InsertNextId(*(polyIds+id));
                vtkIdType list[faceNr];
                vtkIdType replace[faceNr*2];  
                int temp1,temp2;
                //for each face
                for (int pnr = 0; pnr < faceNr; pnr++)
                   {   
                     int lineNr = 0;
                     vtkIdType result = *(faceIds+pnr);
                     geometry->GetPtopeFacets(result,lineIds,lineNr);
                     for ( int q = 0; q < lineNr; q++)
                       { 
                         replace[pnr*2+q] = -1 - *(lineIds+q);                          
                       }
                   }          
               list[0] = replace[0];
               list[1] = replace[1];
               replace[0] = -1;
               replace[1] = -1; 
               // Sort the vertices into right order for a polytope
               for ( int q = 2; q < faceNr; q++)
                  {
                    int value = list[q-1];
                    for ( int p = 1; p < faceNr; p++)
                      {
                        int tempStart = replace[p*2];
                        int tempEnd = replace[p*2+1];
                        if (tempStart > 0)
                          {
                            if (tempStart == value || tempEnd == value ) 
                              { 
                                 if (tempStart == value ) 
                                   { 
                                     list[q] = tempEnd;
                                     replace[p*2]=replace[p*2+1] = -1; 
                                   }
                                 else if (tempEnd == value ) 
                                   {
                                     list[q] = tempStart;
                                     replace[p*2]=replace[p*2+1] = -1;      
                                   }
                              }  
                         }
                    }                     
                } //end vertex order
            
            for (int q = 0; q < faceNr; q++)
                {
                  vtkIdType pid;
                  double p [] = {0,0,0};
                  double rangeP[this->N];
                  geometry->GetActiveDomCoord(-1-list[q], p);
                  geometry->GetActiveRngCoord(-1-list[q], rangeP);
                  // if it is a new point
                  if (pointLocator->InsertUniquePoint(p,pid))
                    {
                     for (int i = 0; i < outFragPD->GetNumberOfArrays(); i++)
                        {
                          outFragPD->GetArray(i)->InsertNextTuple1(rangeP[i]);                               
                        }
                    }
                  poly->InsertNextId(pid);
                }           
             } // End for each polytope  

             outFrags->InsertNextCell(VTK_POLYHEDRON, poly);
             cIdList->InsertNextCell(VTK_POLYGON,polyId);
         }  // End 3D case
       
      } //for each cell
    }  // MAIN LOOP END.

  vtkIdType nrFragments = geometry->GetNrStoredPtopes();
  
  this->UF->SetNumberOfTuples(nrFragments);
  fragSize->SetNumberOfTuples(nrFragments);
  for (vtkIdType i = 0; i < nrFragments; i++)
    {
     this->UF->SetValue(i, -1);
     fragSize->SetValue(i, 1);
    }
  // Compute the quantized slab number for each fragment.
  // Note use of EPSILON to handle problem of numerical approximation
  // pushing slab values to just under a multiple of slab widths, and
  // then use of floor dropping the slab number by one.
  for (int f = 0; f < this->N; f++)
    {
    this->Fields[f]->slabNumbers->SetNumberOfTuples(nrFragments);
    sw   = this->Fields[f]->slabWidth;
    minF = this->Fields[f]->range[0];
    for (vtkIdType i = 0; i < nrFragments; i++)
      {
      fv = geometry->GetStoredRangeComponent(i, f);
      this->Fields[f]->slabNumbers->SetValue(i, floor((fv + EPSILON - minF) / sw));
      }
    }

  // Identify adjacent fragments that should be merged into cells.
  // Iterate over fragment facets, looking for facets that are 
  // shared by two polytopes.  For shared facets, look up the
  // fragment slab numbers.  If these are the same over all
  // components, merge the corresponding union-find partitions.

  vtkIdTypeArray *facets = geometry->GetCenterFacets();
  vtkIdType nrFacets = facets->GetNumberOfTuples();
  
  vtkIdTypeArray *facetsId = geometry->GetCenterFacetsId();
   
  vtkIdType u, v, rootU, rootV;
  bool equal;
  int mergedFragSize;
  // For each fragment facet ...
  for (vtkIdType i = 0; i < nrFacets; i++)
    {
    v = facets->GetComponent(i, 1);
    // Skip if the facet isn't shared.
    if (v < 0) 
      {
      continue;   // no pair of polytopes shares facet i.
      }
    
    u = facets->GetComponent(i, 0);
    
    // u and v are polytopes that share facet i.
    // Check whether the slab number of each component
    // of u matches the corresponding slab number of v.
    equal = true;
    for (int j = 0; equal && j < this->N; j++)
      {
      // Note: slab numbers are short integers: no epsilon involved.
      equal = this->Fields[j]->slabNumbers->GetComponent(u, 0)
              ==
              this->Fields[j]->slabNumbers->GetComponent(v, 0);                
      }
    // If all components are equal, merge the classes.
    if (equal)
      {
      rootU = this->Find(u);
      rootV = this->Find(v);
      if (rootU != rootV)
        {
        this->UF->SetValue(rootU, rootV);
        mergedFragSize = fragSize->GetValue(rootU) + fragSize->GetValue(rootV);
        fragSize->SetValue(rootU, mergedFragSize);
        fragSize->SetValue(rootV, mergedFragSize);
        }
      }
    } // for each fragment facet.
  // Count the number of partitions (-ve entries in the)
  // union-find array.  This is the number of JCN slabs.
  // Keeping a running count of slab number, replace each
  // partition root with an encoding of the final slab id.
  vtkIdType  nrSlabs = 0;
  for (vtkIdType i = 0; i < nrFragments; i++)
    {
    if (this->UF->GetValue(i) < 0)
      {
      // this is the root of a partition.
      this->UF->SetValue(i, -(1 + nrSlabs));      
      nrSlabs++;
      }
    }

 
  // Set up cell vertex table, mapping node ids in full cell
  // graph into node ids in reduced graph.

  cellVertex->SetNumberOfTuples(nrFragments);
  cellVertex->SetName("slab id");
  for (vtkIdType i = 0; i < nrFragments; i++)
    {
      
      rootU = this->UF->GetValue(i);
      while (rootU >= 0)
          {
           rootU = this->UF->GetValue(rootU);
          }
      cellVertex->SetValue(i, -rootU-1);
    }
  // Construct the JCN graph.
  
  // Vertices of the graph correspond to equivalence
  // classes of fragments, i.e. slabs.
    topograph->SetNumberOfVertices(nrSlabs);
    vtkIdType vertU, vertV;
  
  // Two vertices are connected just when there is a face shared by 
  // the two *slabs*, i.e. the slabs are face-adjacent. We again
  // iterate over the facets, but this time for shared facets, we
  // lookup the corresponding slab-id in the union-find structure.
  // If the slab ids are different, the slabs are distinct but
  // adjacent and we create an edge between the corresponding 
  // graph nodes.
  
    
  for (vtkIdType i = 0; i < nrFacets; i++)
    {
    v = facets->GetComponent(i, 1);
    if (v < 0) 
      {
      continue;   // no pair of polytopes shares facet i.
      }
    u = facets->GetComponent(i, 0);

    // Find slab ids for the two fragments
    vertU = this->UF->GetValue(this->Find(u));
    vertV = this->UF->GetValue(this->Find(v));

    // if fragments belong to different slabs, add an
    // edge between the nodes representing the polytopes.
    // Note the conversion from encoded slab number in the
    // UF structure back to slab ids.

    if (vertU != vertV)
      {
      u = -vertU - 1;
      v = -vertV - 1;
      if (topograph->GetEdgeId(u, v) < 0)
        {
        topograph->AddEdge(u, v);
        }
      }
    }  // for each shared facet / potential edge.
  vtkIdType topoU, topoV,size;
  double p[MaxDomainDim];  
  p[0] = p[1] = p[2] = 0; 

  // Identify boundary slabs.  Brute force, boundary flags are initialised 
  // to false.  We scan facets looking for a facet with only one neighbouring 
  // polytope.  The boundary flag for this polytope is set to true - one
  // boundary facet is sufficient.

  // Variable holding array is declared at function level as we need to
  // refer to it later at the end of the method for setting up output.
  vtkSmartPointer<vtkBitArray> isOnBoundary 
      = vtkSmartPointer<vtkBitArray>::New();

  if (this->IdentifyBoundarySlabs)
    {
    // Create array for holding boundary flag.
    isOnBoundary->SetName("Boundary");
    isOnBoundary->SetNumberOfValues(nrSlabs);
    // Initialise entries to false - assume all are internal.
    for (vtkIdType i = 0; i < nrSlabs; i++)
      {
      isOnBoundary->SetValue(i, 0);
      }
    // Check each facet.
    for (vtkIdType i = 0; i < nrFacets; i++)
      {
      v = facets->GetComponent(i, 1);
      // Is there a second entry for the facet?
      if (v >= 0) 
        {
        continue;   // part of a pair, so not a boundary facet.
        }
      // If only one entry, that is the fragment ID of the polytope
      // that lies on the boundary.  Use the union-find structure 
      // to find the corresponding node id.
      u = facets->GetComponent(i, 0);
      vertU = this->UF->GetValue(this->Find(u));
      // Convert from root marker to node id.
      vertU = -vertU - 1;
      isOnBoundary->SetValue(vertU, 1);
      }
    }
    

  // Copy point data from the equivalence class of the slab into entry for 
  // the corresponding graph vertex.
  // It is more efficient to do this as a separate step than merge into
  // an earlier loop, as we can now pre-allocate space in the node array 
  // and avoid resizing.
  
  for (int j = 0; j < this->N; j++)
    {
    this->Fields[j]->nodeScalars->SetNumberOfTuples(nrSlabs);
    }  
    
  // Initialise node size array.
  nodeSize->SetNumberOfTuples(nrSlabs);
  nodeSize->SetName("size");
  // NOTE - may still be more efficient to swap the loops and load
  // slab width and base once into variables, rather than looking up
  // in the field.
  for (vtkIdType i = 0; i < nrFragments; i++)
    {
    vtkIdType s = this->UF->GetValue(i);
    if (s < 0)
      {
      // decode slab id
      s = -s - 1;
      // compute the slab value of each field from the slab number.
      for (int j = 0; j < this->N; j++)
        {
        this->Fields[j]->nodeScalars->SetComponent(s, 0,
            this->Fields[j]->range[0] 
            + 
            this->Fields[j]->slabWidth * this->Fields[j]->slabNumbers->GetComponent(i, 0));
        }
       nodeSize->SetValue(s, fragSize->GetValue(i)+1); 
      }
    } // Copy global scalars to node scalars.

   // Copy the JCN graph to the output.
   output->ShallowCopy(topograph);    
   outFrags->SetPoints(newPoints);
   outFrags->BuildLinks();
  for (std::vector<Field*>::iterator it = this->Fields.begin();
      it != this->Fields.end();
      it++)
    {
    outFragCD->AddArray((*it)->slabNumbers);
    }
  
 
  outFragCD->AddArray(cellVertex);
  outFrags->Squeeze();    

   // boundary computation
  vtkIdType faceNr;
  int polyNr;
  vtkIdType bid;
  vtkIdType cellU, cellV;
  vtkIdType thisNodeIx, nextNodeIx;
  vtkSmartPointer<vtkIdList> line  = vtkSmartPointer<vtkIdList>::New();
  line->SetNumberOfIds(2);
  vtkIdType *polyIds = NULL;
  
  vtkSmartPointer<vtkIdList> lineIds  = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cell  = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> facePtr  = vtkSmartPointer<vtkIdList>::New();
   // Global stats
  vtkIdType bordTriCount = 0;
  vtkIdType pId[2];

  vtkIdList *(*borders) = new vtkIdList*[nrSlabs];
  for (vtkIdType i = 0; i < nrSlabs; i++)
    {
    borders[i] = vtkIdList::New();
    if (this->M == 3)
      {
      // Add counter for number of polygons.
      borders[i]->InsertNextId(0);
      }
    }

  // Compute slab geometry
  // facets : list of edges
  if (this->M == 2)
    {
    // Number of edges
    for (vtkIdType i = 0; i < nrFacets; i++)
       {
       // u: which fragment the edge i is in. 
       u = facets->GetComponent(i,0);
       // v: which fragment the edge i's neighbour is in. 
       v = facets->GetComponent(i,1);
       // which slab it is in
       cellU = cellVertex->GetValue(u);
           
       // Get the edge end points' ids
       pId[0] = facetsId->GetComponent(i, 0);
       pId[1] = facetsId->GetComponent(i,1);  

       // If the edge is on the boundary
       if (v < 0)
         {        
         // Outer edge of the dataset, so automatically
         // included in the region boundary.
         outFrags->GetCellPoints(u, line);
         cIdList->GetCellPoints(u, lineIds);
         for (int pnr = 0; pnr < 2; pnr++)
            {   
            int index =  pId[pnr];
            for (int m = 0; m < lineIds->GetNumberOfIds(); m++)
               {
               if (lineIds->GetId(m) == index) 
                 {
                  borders[cellU]->InsertNextId(line->GetId(m));
                  break;
                  }   
                }                         
             } 
           bordTriCount += 1;
           } 
       else 
         {
         cellV = cellVertex->GetValue(v);
         if (cellU == cellV)
           {
           continue;
           }
         // edge shared by two slabs, then add this edge to each of the shared slab
         outFrags->GetCellPoints(u, line);
         cIdList->GetCellPoints(u, lineIds);
              
         for (int pnr = 0; pnr < 2; pnr++)
            {   
            int index =  pId[pnr];
            for (int m = 0; m < lineIds->GetNumberOfIds(); m++)
               {
               if (lineIds->GetId(m) == index) 
                 {
                  borders[cellU]->InsertNextId(line->GetId(m));
                  break;
                 }   
                }                         
              } 

         pId[0] = facetsId->GetComponent(i, 2);
         pId[1] = facetsId->GetComponent(i,3);    
         outFrags->GetCellPoints(v, line);
         cIdList->GetCellPoints(v, lineIds);
              
         for (int pnr = 0; pnr < 2; pnr++)
            {   
            int index =  pId[pnr];
            for (int m = 0; m < lineIds->GetNumberOfIds(); m++)
               {
               if (lineIds->GetId(m) == index) 
                 {
                 borders[cellV]->InsertNextId(line->GetId(m));
                 break;
                 }   
               }                         
            }
          bordTriCount += 2;
          }
        } 
       
        // Convert list of vertex pairs into a
        // single polygon.  Painful.  We make an 
        // inverse table to track each of the two
       // indices in the list where a given vertex
       // appears, and then play joint-the-dots.
       for (vtkIdType i = 0; i < nrSlabs; i++)
          { 
          vtkIdList *poly = borders[i];         
          vtkIdType baseId, lastId;
          baseId = lastId = poly->GetId(0);
          for (vtkIdType j = 0; j < poly->GetNumberOfIds(); j++)
             {
             vtkIdType id = poly->GetId(j);
             if (id < baseId)
               {
                baseId = id;
               }
              if (id > lastId)
               {
                lastId = id;
               }
             }

          const vtkIdType invSize = 2 * (lastId-baseId+1);
          vtkIdType *inv    = new vtkIdType [invSize];

          for (vtkIdType j = 0; j < invSize; j++)
             {
             inv[j] = -1;
             }
            
          for (vtkIdType j = 0; j < poly->GetNumberOfIds(); j++)
             {
             vtkIdType offset = 2*(poly->GetId(j) - baseId);

             if (inv[offset] < 0)
               {
               inv[offset] = j;
               }
             else
               {
               inv[offset + 1] = j;
               }
              }

           // Now join-the-dots.
           line->Reset();
           thisNodeIx = 0;
           do {      
              bid = poly->GetId(thisNodeIx);       
              line->InsertNextId(bid);          
              nextNodeIx = inv[2*(bid - baseId)];
              if (nextNodeIx == thisNodeIx)
                {
                 nextNodeIx = inv[2*(bid - baseId)+1];
                }
              // Select the "other" node on the edge.
              thisNodeIx = (nextNodeIx % 2) ? nextNodeIx-1 : nextNodeIx+1;
             } // do  
           while (thisNodeIx);
	
           // ... and insert polygon into output
           outBoundary->InsertNextCell(VTK_POLYGON, line);
           delete [] inv;
      } // For each vertex 
    } // Dimension 2
  else if (M==3)
     {
      // Number of faces
      for (vtkIdType i = 0; i < nrFacets; i++)
         { 
         // u: which fragment the facet i is in. 
         u = facets->GetComponent(i,0);
         v = facets->GetComponent(i,1);
         cellU = cellVertex->GetValue(u);
         bid = facetsId->GetComponent(i,0);
         cIdList->GetCellPoints(u,line);
         // Get the id of the current facet. 
         for (int m = 0; m < line->GetNumberOfIds(); m++)
            {
            if (line->GetId(m) == bid)
              {
              bid = m;
              break;
              }
             }   
          if (v < 0)
            {
            // Outer edge of the dataset, so automatically
            // included in the region boundary.
            facePtr = outFrags->GetCell(u)->GetFace(bid)->GetPointIds();
            cell = borders[cellU];
            cell->InsertNextId(facePtr->GetNumberOfIds());
            for (int j = 0; j < facePtr->GetNumberOfIds(); j++)
               {
               cell->InsertNextId(facePtr->GetId(j));
               }
            // increment polygon count
            cell->SetId(0, cell->GetId(0) + 1);
            } // if outer edge (v < 0)       
          else
            {
            // Check whether the two cells belong to the same 
            // region.  If the same, this is an interior border
            // and we skip over it.
            cellV = cellVertex->GetValue(v);
            if (cellU == cellV)
              {
              continue;
              }
             // u and v are different regions, so add the
             // boundary to both region polyhedra.
             facePtr = outFrags->GetCell(u)->GetFace(bid)->GetPointIds();
             cell = borders[cellU];
             cell->InsertNextId(facePtr->GetNumberOfIds());
             for (int j = 0; j < facePtr->GetNumberOfIds(); j++)
                {
                cell->InsertNextId(facePtr->GetId(j));
                }
             // increment polygon count for cell U
            cell->SetId(0, cell->GetId(0) + 1);

            bid = facetsId->GetComponent(i,1);
            cIdList->GetCellPoints(v,line);
            // Obtain the current facet index
            for (int m = 0; m < line->GetNumberOfIds(); m++)
               {
               if (line->GetId(m) == bid)
                 {
                 bid = m;
                 break;
                 }
               }
                   
            facePtr = outFrags->GetCell(v)->GetFace(bid)->GetPointIds();
            cell = borders[cellV];
            cell->InsertNextId(facePtr->GetNumberOfIds());
            for (int j = 0; j < facePtr->GetNumberOfIds(); j++)
               {
               cell->InsertNextId(facePtr->GetId(j));
               }
             //increment polygon count for cell V
             cell->SetId(0, cell->GetId(0) + 1);
             }  // shared inner edge
       }
      
      for (vtkIdType i = 0; i < nrSlabs; i++)
         {
          outBoundary->InsertNextCell(VTK_POLYHEDRON, borders[i]);  
          bordTriCount += borders[i]->GetId(0);
         }    
  }
  outBoundary->SetPoints(newPoints);
  outBoundary->BuildLinks();
  for (vtkIdType i = 0; i < nrSlabs; i++)
     {
     borders[i]->Delete();
     }
  delete [] borders;
  outBoundary->Squeeze();
   
  // Initialize union-find structure to one entry per fragment, 
  // using -1 to indicate partition roots.
  const double origin[] = {0.0, 0.0, 0.0};
  // Initialise array for tracking fragment locations
  fragCenters->SetNumberOfComponents(3);
  fragCenters->SetNumberOfTuples(nrSlabs);
  fragCenters->SetName("center"); 
  slabSize->SetNumberOfTuples(nrSlabs);
   
  for (vtkIdType i = 0; i < nrSlabs; i++)
     {
     fragCenters->SetTupleValue(i,origin);
      slabSize->SetValue(i,0);
     }
    
  // Compute the center point of each fragment
  for (vtkIdType i = 0; i < nrFacets; i++)
     {
      geometry->GetCenterLocator()->GetPoint(i,p);
    
      v = facets->GetComponent(i, 1);
      u = facets->GetComponent(i, 0);
     
      topoU = cellVertex->GetValue(u);
      slabSize->SetValue(topoU, slabSize->GetValue(topoU)+1);
      for (int j = 0; j < 3; j++)
         {
         fragCenters->SetComponent(topoU, j, fragCenters->GetComponent(topoU, j)+p[j]);
         }  
       // Skip if the facet isn't shared.
      if (v < 0) 
        {
        continue;   // no pair of polytopes shares facet i.
        }  
      topoV = cellVertex->GetValue(v); 
      if (topograph->GetEdgeId(topoU, topoV) >= 0)
        {
        continue;
        }
      for (int j = 0; j < 3; j++)
        {
        fragCenters->SetComponent(topoV, j, fragCenters->GetComponent(topoV, j)+p[j]);
        } 
       slabSize->SetValue(topoV, slabSize->GetValue(topoV)+1);
      } // for each fragment facet.
    

  // Complete computation of slab centers.
  for (vtkIdType u = 0; u < nrSlabs; u++)
    {
     size = slabSize->GetValue(u);
     for (int j = 0; j < 3; j++)
       {
       fragCenters->SetComponent(u, j, fragCenters->GetComponent(u, j) / size);
       } 
    }
    
  // Output points (vertex locations) have not been defined, but some downstream filters
  // get upset if no point object is present, so we simply copy the input points across,
  // in the same form that they arrived.
  if (this->UseSpatialCoordinateArray)
    {
    output->GetVertexData()->AddArray(points);
    }
  else
    {
    output->SetPoints(input->GetPoints());
    }

  // Add vertex data arrays to the output, and deallocate working arrays.
  for (int j = 0; j < this->N; j++)
    {
    output->GetVertexData()->AddArray(this->Fields[j]->nodeScalars);
    this->Fields[j]->slabNumbers = NULL;
    this->Fields[j]->nodeScalars = NULL;
    }
  if (this->IdentifyBoundarySlabs)
    {
    output->GetVertexData()->AddArray(isOnBoundary);
    }
  output->GetVertexData()->AddArray(nodeSize);
  output->GetVertexData()->AddArray(fragCenters);
  return 1; 
}

