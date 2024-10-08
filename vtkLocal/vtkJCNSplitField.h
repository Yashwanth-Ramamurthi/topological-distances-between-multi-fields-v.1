/*=========================================================================

 *	File: vtkJCNSplitField.h
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
// .NAME vtkJCNSplitField - Split a field into single component fields
// .SECTION Description
// vtkJCNSplitField is used to split a multi-component field (vtkDataArray)
// into multiple single component fields. The new fields are put in
// the same field data as the original field. The output arrays
// are of the same type as the input array. Example:
// @verbatim
// sf->SetInputField("gradient", vtkJCNSplitField::POINT_DATA);
// sf->Split(0, "firstcomponent");
// @endverbatim
// tells vtkJCNSplitField to extract the first component of the field
// called gradient and create an array called firstcomponent (the
// new field will be in the output's point data).
// The same can be done from Tcl:
// @verbatim
// sf SetInputField gradient POINT_DATA
// sf Split 0 firstcomponent
//
// AttributeTypes: SCALARS, VECTORS, NORMALS, TCOORDS, TENSORS
// Field locations: DATA_OBJECT, POINT_DATA, CELL_DATA
// @endverbatim
// Note that, by default, the original array is also passed through.

// .SECTION Caveats
// When using Tcl, Java, Python or Visual Basic bindings, the array name 
// can not be one of the  AttributeTypes when calling Split() which takes
// strings as arguments. The Tcl (Java etc.) command will
// always assume the string corresponds to an attribute type when
// the argument is one of the AttributeTypes. In this situation,
// use the Split() which takes enums.

// .SECTION See Also
// vtkFieldData vtkDataSet vtkDataObjectToDataSetFilter
// vtkDataSetAttributes vtkDataArray vtkRearrangeFields
// vtkAssignAttribute vtkJCNMergeFields

#ifndef __vtkJCNSplitField_h
#define __vtkJCNSplitField_h

#include "vtkPassInputTypeAlgorithm.h"

#include "vtkDataSetAttributes.h" // Needed for NUM_ATTRIBUTES
#include "vtkFiltersGeneralModule.h" // For export macro


class vtkDataArray;
class vtkFieldData;

class VTKFILTERSGENERAL_EXPORT vtkJCNSplitField : public vtkPassInputTypeAlgorithm
{
public:
  vtkTypeMacro(vtkJCNSplitField,vtkPassInputTypeAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Create a new vtkJCNSplitField.
  static vtkJCNSplitField *New();

  // Description:
  // Use the  given attribute in the field data given
  // by fieldLoc as input.
  void SetInputField(int attributeType, int fieldLoc);

  // Description:
  // Use the array with given name in the field data given
  // by fieldLoc as input.
  void SetInputField(const char* name, int fieldLoc);

  // Description:
  // Helper method used by other language bindings. Allows the caller to
  // specify arguments as strings instead of enums.
  void SetInputField(const char* name, const char* fieldLoc);

  // Description:
  // Create a new array with the given component.
  void Split(int component, const char* arrayName);

//BTX
  enum FieldLocations
  {
    DATA_OBJECT=0,
    POINT_DATA=1,
    CELL_DATA=2,
    VERTEX_DATA=3,
    EDGE_DATA=4
  };
//ETX

//BTX
  struct Component
  {
    int Index;
    char* FieldName;   
    Component* Next;   // linked list
    void SetName(const char* name)
      {
        delete[] this->FieldName;
        this->FieldName = 0;
        if (name)
          {
          this->FieldName = new char[strlen(name)+1];
          strcpy(this->FieldName, name);
          }
      }
    Component() { FieldName = 0; }
    ~Component() { delete[] FieldName; }
  };
//ETX

protected:

//BTX
  enum FieldTypes
  {
    NAME,
    ATTRIBUTE
  };
//ETX

  vtkJCNSplitField();
  virtual ~vtkJCNSplitField();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  char* FieldName;
  int FieldType;
  int AttributeType;
  int FieldLocation;

  static char FieldLocationNames[5][12];
  static char AttributeNames[vtkDataSetAttributes::NUM_ATTRIBUTES][10];

  vtkDataArray* SplitArray(vtkDataArray* da, int component);


  // Components are stored as a linked list.
  Component* Head;
  Component* Tail;

  // Methods to browse/modify the linked list.
  Component* GetNextComponent(Component* op)
    { return op->Next; }
  Component* GetFirst()
    { return this->Head; }
  void AddComponent(Component* op);
  Component* FindComponent(int index);
  void DeleteAllComponents();

  void PrintComponent(Component* op, ostream& os, vtkIndent indent);
  void PrintAllComponents(ostream& os, vtkIndent indent);
private:
  vtkJCNSplitField(const vtkJCNSplitField&);  // Not implemented.
  void operator=(const vtkJCNSplitField&);  // Not implemented.
};

#endif


