/*=========================================================================

 *	File: vtkJCNMergePoints.h
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
// .NAME vtkJCNMergePoints - merge exactly coincident points
// .SECTION Description
// vtkJCNMergePoints is a locator object to quickly locate points in 3D.
// The primary difference between vtkJCNMergePoints and its superclass
// vtkPointLocator is that vtkJCNMergePoints merges precisely coincident points
// and is therefore much faster.
// .SECTION See Also
// vtkCleanPolyData

#ifndef __vtkJCNMergePoints_h
#define __vtkJCNMergePoints_h

#include "vtkCommonDataModelModule.h" // For export macro
#include "vtkPointLocator.h"

class VTKCOMMONDATAMODEL_EXPORT vtkJCNMergePoints : public vtkPointLocator
{
public:
  static vtkJCNMergePoints *New();
  vtkTypeMacro(vtkJCNMergePoints,vtkPointLocator);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Determine whether point given by x[3] has been inserted into points list.
  // Return id of previously inserted point if this is true, otherwise return
  // -1.
  vtkIdType IsInsertedPoint(const double x[3]);
  vtkIdType IsInsertedPoint(double x, double  y, double z)
    {return this->vtkPointLocator::IsInsertedPoint(x, y, z); };

  // Description:
  // Determine whether point given by x[3] has been inserted into points list.
  // Return 0 if point was already in the list, otherwise return 1. If the
  // point was not in the list, it will be ADDED.  In either case, the id of
  // the point (newly inserted or not) is returned in the ptId argument.
  // Note this combines the functionality of IsInsertedPoint() followed
  // by a call to InsertNextPoint().
  int InsertUniquePoint(const double x[3], vtkIdType &ptId);
  
protected:
  vtkJCNMergePoints() {};
  ~vtkJCNMergePoints() {};
  
private:
  vtkJCNMergePoints(const vtkJCNMergePoints&);  // Not implemented.
  void operator=(const vtkJCNMergePoints&);  // Not implemented.
};

#endif


