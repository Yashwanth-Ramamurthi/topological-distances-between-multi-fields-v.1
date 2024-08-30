/*=========================================================================
This source has no copyright.  It is intended to be copied by users
wishing to create their own VTK classes locally.
=========================================================================*/

// .NAME vtkPyramidTree - Data structure for searching in high-dimensional space.
// .SECTION Description
// This class implements a Pyramid Technique.  The Pyramid-tree is based on a
// partitioning strategy which is optimized for high-dimensional data. The basic 
// idea is to divide the data space first into 2d pyramids sharing the
// center point of the space as a top.  In a second step, the
// single pyramids are cut into slices parallel to the basis
// of the pyramid. These slices form the data pages. Furthermore,
// this partition provides a mapping from the given d-dimensional space to a 
// 1-dimensional space. We will use hashmap on this 1-dimensional space.


#ifndef __vtkPyramidTree_h
#define __vtkPyramidTree_h

#include "vtkMETAModule.h"
#include "vtkObject.h"
#include <vector>
#include <boost/unordered_map.hpp>

// Description:
// Small value to test the equality between two values.

#define MYEPSILON 1.0e-7
#define MAX_BUCKET_NUMBER 300000


// Equality measure for hash map function
struct iequal_to : std::binary_function<double, double, bool>
{
  bool operator()(double const& x, double const& y) const
    {
    if (fabs(x - y) < MYEPSILON)
      {
      return 1;
      }
     else
      {
       return 0;
      }
    }
};

// hash function used in hash map
struct ihash
: std::unary_function<int, size_t>
{
  size_t operator()(int const& x) const
    {
    return x;    
    }
};


class VTK_META_EXPORT vtkPyramidTree : public vtkObject
{
public:
  static vtkPyramidTree* New();
  vtkTypeMacro(vtkPyramidTree, vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);
  
  // Description: 
  // Set the dimensionality of the data set
  void SetDimension(int dim);

  // Description:  
  // Set the minimal bounds of each dimension
  void SetMinBounds(const double* minBounds);

  // Description:  
  // Set the maximal bounds of each dimension
  void SetMaxBounds(const double* maxBounds);

  // Description:  
  // Return the number of elements of the inserted data.
  int GetNumberElement();

  // Description:  
  // Return the number of dimensions of the inserted data. 
  int GetNumberDimension();

  // Description: 
  // Insert an (m-dimensional) point into the structure.  If point is already present in
  //  the structure, return 0, and set pid to that point's unique id.  If the point was
  //  not in the structure, return 1 and set the pid to the id allocated to the point.
  //  The library can assume that "point" points to at least dim doubles, where dim was
  //  the value provided in SetDimension.
  int InsertUniquePoint(const double *point, vtkIdType &pid);

  // Description: 
  // Retrieve a point with a given id
  void GetPointIndex(vtkIdType &pid, const double *point);

  // Description:  
  // Initialize the data structure
  void Initialize(int col, const double *bounds);

  // Description:  
  // Given a point ID, retrive the actual points
  void GetPoint(vtkIdType pid, double *point);

  // Description:
  // Return the number of points in the data structure
  int GetNumberOfPoints();

  // Description:  
  // Insert one point into the data structure. 
  void InsertNextPoint(const double *p);

  // Description:
  // Time for performance testing.
  double GetTime();

  // Description:
  // Return number of bucket.
  int GetBucketCount();

  // Description:
  // Return the number of elements.
  int GetBucketSize();
 
protected:
  vtkPyramidTree();
  ~vtkPyramidTree();

private:  
  // Insert a data point into data strucuture
  void InsertToStructure(std::vector<double> myPoint, vtkIdType dupliate );

  // Compute PT Structure, return the key value
  int ComputePTStruture(std::vector<double> myPoint);

  // Hash function
  int HashValue(std::vector<double> myPoint);

  //Original Data Set
  std::vector< std::vector<double> > RawData;

  // Dimension and size of the data set.
  int ColNumber;
 
  // sum of the element in each row.  For testing the overlapping keys in the future.
  std::vector<double> SumData;

  // Maximal bounds of the data set.
  std::vector<double> MaxDim;

  // Mimial bounds of the data set.
  std::vector<double> MinDim;
 
  // Median values of the data set.
  std::vector<double> Median;
 
  // Unordered_map for storing the key values of 1-dimensional interpolated point. 
  boost::unordered_map<int, std::vector<int>,ihash> Map;

  // Time for performance testing.
  double TTime;

  // Bucket Interval.
  int Interval;
};

// Helper function:
// Compare two double values subject to an error tolerance.
// Return -1, 0, 1 respectively if t is less than,
// approximately equal to, or greater than, base.
inline int compare(double t, double base)
{
  if (fabs(t - base) < MYEPSILON)
    {
    return 0;
    }
  if (t < base)
    {
    return -1;
    }
  return 1;
}

// Compare two vectors according to the difference of each dimension of them.
inline int compareVector(std::vector<double> a, std::vector<double> b)
{
  int sizeA = a.size();
  int sizeB = b.size();
    
  if (sizeA == sizeB)
    {
    for (int i = 0; i < sizeA; i++)
       {
       // if two values are not equal
       if (compare(a.at(i),b.at(i)))
         {
         return 1;
         break;
         }
       } 
    }
  else
    {
    //not equal return 1
    return 1;
    }
  //equal return 0
  return 0;
}

#endif
