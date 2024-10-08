/*=========================================================================

 *	File: vtkJCNMergeFields.h
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
// .NAME vtkJCNMergeFields - Merge multiple fields into one.
// .SECTION Description
// vtkJCNMergeFields is used to merge mutliple field into one.
// The new field is put in the same field data as the original field.
// For example
// @verbatim
// mf->SetOutputField("foo", vtkJCNMergeFields::POINT_DATA);
// mf->SetNumberOfComponents(2);
// mf->Merge(0, "array1", 1);
// mf->Merge(1, "array2", 0);
// @endverbatim
// will tell vtkJCNMergeFields to use the 2nd component of array1 and
// the 1st component of array2 to create a 2 component field called foo.
// The same can be done using Tcl:
// @verbatim
// mf SetOutputField foo POINT_DATA
// mf Merge 0 array1 1
// mf Merge 1 array2 0
//
// Field locations: DATA_OBJECT, POINT_DATA, CELL_DATA
// @endverbatim

// .SECTION See Also
// vtkFieldData vtkDataSet vtkDataObjectToDataSetFilter
// vtkDataSetAttributes vtkDataArray vtkRearrangeFields
// vtkSplitField vtkAssignAttribute

#ifndef __vtkJCNMergeFields_h
#define __vtkJCNMergeFields_h

#include "vtkPassInputTypeAlgorithm.h"
#include "vtkFiltersCoreModule.h" // For export macro

class vtkDataArray;
class vtkFieldData;

class VTKFILTERSCORE_EXPORT vtkJCNMergeFields : public vtkPassInputTypeAlgorithm
{
public:
  vtkTypeMacro(vtkJCNMergeFields,vtkPassInputTypeAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Create a new vtkJCNMergeFields.
  static vtkJCNMergeFields *New();

  // Description:
  // The output field will have the given name and it will be in
  // fieldLoc (the input fields also have to be in fieldLoc).
  void SetOutputField(const char* name, int fieldLoc);

  // Description:
  // Helper method used by the other language bindings. Allows the caller to
  // specify arguments as strings instead of enums.Returns an operation id 
  // which can later be used to remove the operation.
  void SetOutputField(const char* name, const char* fieldLoc);

  // Description:
  // Add a component (arrayName,sourceComp) to the output field.
  void Merge(int component, const char* arrayName, int sourceComp);

  // Description:
  // Set the number of the components in the output field.
  // This has to be set before execution. Default value is 0.
  vtkSetMacro(NumberOfComponents, int);
  vtkGetMacro(NumberOfComponents, int);

  enum FieldLocations
  {
    DATA_OBJECT=0,
    POINT_DATA=1,
    CELL_DATA=2,
    VERTEX_DATA=3,
    EDGE_DATA=4
  };

//BTX
  struct Component
  {
    int Index;
    int SourceIndex;
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
  enum FieldType
  {
    NAME,
    ATTRIBUTE
  };
//ETX

  vtkJCNMergeFields();
  virtual ~vtkJCNMergeFields();

  int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

  char* FieldName;
  int FieldLocation;
  int NumberOfComponents;
  int OutputDataType;

  static char FieldLocationNames[5][12];


  int MergeArray(vtkDataArray* in, vtkDataArray* out, int inComp, int outComp);

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
  vtkJCNMergeFields(const vtkJCNMergeFields&);  // Not implemented.
  void operator=(const vtkJCNMergeFields&);  // Not implemented.
};

#endif


