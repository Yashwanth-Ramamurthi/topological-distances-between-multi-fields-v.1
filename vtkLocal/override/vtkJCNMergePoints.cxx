/*=========================================================================

 *	File: vtkJCNMergePoints.cxx
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
#include "vtkJCNMergePoints.h"

#include "vtkDataArray.h"
#include "vtkIdList.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkFloatArray.h"

vtkStandardNewMacro(vtkJCNMergePoints);

// Determine whether point given by x[3] has been inserted into points list.
// Return id of previously inserted point if this is true, otherwise return
// -1.
vtkIdType vtkJCNMergePoints::IsInsertedPoint(const double x[3])
{
  vtkIdType i, ijk0, ijk1, ijk2;
  vtkIdType idx;
  vtkIdList *bucket;
  //
  //  Locate bucket that point is in.
  //
  ijk0 = static_cast<vtkIdType>(
    static_cast<double> ((x[0] - this->Bounds[0]) / 
                         (this->Bounds[1] - this->Bounds[0]))
    * (this->Divisions[0] - 1));
  ijk1 = static_cast<vtkIdType>(
    static_cast<double> ((x[1] - this->Bounds[2]) / 
                         (this->Bounds[3] - this->Bounds[2]))
    * (this->Divisions[1] - 1));
  ijk2 = static_cast<vtkIdType>(
    static_cast<double> ((x[2] - this->Bounds[4]) / 
                         (this->Bounds[5] - this->Bounds[4]))
    * (this->Divisions[2] - 1));


  idx = ijk0 + ijk1*this->Divisions[0] + 
        ijk2*this->Divisions[0]*this->Divisions[1];

  bucket = this->HashTable[idx];

  if ( ! bucket )
    {
    return -1;
    }
  else // see whether we've got duplicate point
    {
    //
    // Check the list of points in that bucket.
    //
    vtkIdType ptId;
    int nbOfIds = bucket->GetNumberOfIds ();

    // For efficiency reasons, we break the data abstraction for points
    // and ids (we are assuming and vtkIdList
    // is storing ints).
    vtkDataArray *dataArray = this->Points->GetData();
    vtkIdType *idArray = bucket->GetPointer(0);
    if (dataArray->GetDataType() == VTK_FLOAT)
      {
      float f[3];
      f[0] = static_cast<float>(x[0]);
      f[1] = static_cast<float>(x[1]);
      f[2] = static_cast<float>(x[2]);
      vtkFloatArray *floatArray = static_cast<vtkFloatArray *>(dataArray);
      float *pt;
      for (i=0; i < nbOfIds; i++) 
        {
        ptId = idArray[i];
        pt = floatArray->GetPointer(0) + 3*ptId;
        if ( f[0] == pt[0] && f[1] == pt[1] && f[2] == pt[2] )
          {
          return ptId;
          }
        }
      }
    else
      {
      // Using the double interface
      double *pt;
      for (i=0; i < nbOfIds; i++) 
        {
        ptId = idArray[i];
        pt = dataArray->GetTuple(ptId);
        if ( x[0] == pt[0] && x[1] == pt[1] && x[2] == pt[2] )
          {
          return ptId;
          }
        }
      }
    }

  return -1;
}


int vtkJCNMergePoints::InsertUniquePoint(const double x[3], vtkIdType &id)
{
  vtkIdType i, ijk0, ijk1, ijk2;
  vtkIdType idx;
  vtkIdList *bucket;

  //
  //  Locate bucket that point is in.
  //
  ijk0 = static_cast<vtkIdType>(
    static_cast<double> ((x[0] - this->Bounds[0]) / 
                         (this->Bounds[1] - this->Bounds[0]))
    * (this->Divisions[0] - 1));
  ijk1 = static_cast<vtkIdType>(
    static_cast<double> ((x[1] - this->Bounds[2]) / 
                         (this->Bounds[3] - this->Bounds[2]))
    * (this->Divisions[1] - 1));
  ijk2 = static_cast<vtkIdType>(
    static_cast<double> ((x[2] - this->Bounds[4]) / 
                         (this->Bounds[5] - this->Bounds[4]))
    * (this->Divisions[2] - 1));

  idx = ijk0 + ijk1*this->Divisions[0] + 
        ijk2*this->Divisions[0]*this->Divisions[1];

/*
cout << endl << "@@ " << x[0] << "," << x[1] << "," << x[2];
cout << "      " << ijk0 << " " << ijk1 << " " << ijk2 << " " << idx << endl;
*/

  bucket = this->HashTable[idx];

  if (bucket) // see whether we've got duplicate point
    {
    //
    // Check the list of points in that bucket.
    //
    vtkIdType ptId;
    int nbOfIds = bucket->GetNumberOfIds ();

    // For efficiency reasons, we break the data abstraction for points
    // and ids (we are assuming vtkPoints stores a vtkIdList
    // is storing ints).
    vtkDataArray *dataArray = this->Points->GetData();
    vtkIdType *idArray = bucket->GetPointer(0);
    
    if (dataArray->GetDataType() == VTK_FLOAT)
      {
      float f[3];
      f[0] = static_cast<float>(x[0]);
      f[1] = static_cast<float>(x[1]);
      f[2] = static_cast<float>(x[2]);
      vtkFloatArray *floatArray = static_cast<vtkFloatArray *>(dataArray);
      float *pt;
      unsigned *ptf = (unsigned*)pt;
      for (i=0; i < nbOfIds; i++) 
        {
        ptId = idArray[i];
        pt = floatArray->GetPointer(0) + 3*ptId;
//        if ( f[0] == pt[0] && f[1] == pt[1] && f[2] == pt[2] )
/*
cout << "$$F " << f[0] << "," << f[1] << "," << f[2];
cout << "   " << pt[0] << "," << pt[1] << "," << pt[2];
cout << "   " << (f[0] - pt[0]) << " " << (f[1] - pt[1]) << " " << (f[2] - pt[2]);
*/
        if ( fabsf(f[0] - pt[0]) < 0.000001 
          && fabsf(f[1] - pt[1]) < 0.000001 
          && fabsf(f[2] - pt[2]) < 0.000001 )
/*
        if ( fabsf(f[0] - pt[0]) < 0.00000001 
          && fabsf(f[1] - pt[1]) < 0.00000001 
          && fabsf(f[2] - pt[2]) < 0.00000001 )
*/
          {
          // point is already in the list, return 0 and set the id parameter
//cout << "OK " << ptId << endl;
          id = ptId;
          return 0;
          }
//cout << "FAIL" << endl;          
        }
      }
    else
      {
      // Using the double interface
      double *pt;
      for (i=0; i < nbOfIds; i++) 
        {
        ptId = idArray[i];
        pt = dataArray->GetTuple(ptId);

/*
cout << "$$D " << x[0] << "," << x[1] << "," << x[2];
cout << "   " << pt[0] << "," << pt[1] << "," << pt[2];
cout << "   " << (x[0] - pt[0]) << " " << (x[1] - pt[1]) << " " << (x[2] - pt[2]);
*/
/*
        if ( fabs(x[0] - pt[0]) < 0.000001 
          && fabs(x[1] - pt[1]) < 0.000001 
          && fabs(x[2] - pt[2]) < 0.000001 )
*/
        if ( fabs(x[0] - pt[0]) < 0.0000001 
          && fabs(x[1] - pt[1]) < 0.0000001 
          && fabs(x[2] - pt[2]) < 0.0000001 )
//        if ( x[0] == pt[0] && x[1] == pt[1] && x[2] == pt[2] )
          {
          // point is already in the list, return 0 and set the id parameter
          id = ptId;
//cout << "OK " << ptId << endl;
          return 0;
          }
//cout << "FAIL" << endl;          
          
        }
      }
    }
  else
    {
    // create a bucket point list and insert the point
    bucket = vtkIdList::New();
    bucket->Allocate(this->NumberOfPointsPerBucket/2,
                     this->NumberOfPointsPerBucket/3);
    this->HashTable[idx] = bucket;
    }

  // point has to be added
  bucket->InsertNextId(this->InsertionPointId);
  this->Points->InsertPoint(this->InsertionPointId,x);
  id = this->InsertionPointId++;

  return 1;
}

//----------------------------------------------------------------------------
void vtkJCNMergePoints::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}
