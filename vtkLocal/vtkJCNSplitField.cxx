/*=========================================================================

 *	File: vtkJCNSplitField.cxx
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/
#include "vtkJCNSplitField.h"

#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkDataObject.h"
#include "vtkGraph.h"
#include "vtkDataSetAttributes.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include <ctype.h>

vtkStandardNewMacro(vtkJCNSplitField);

char vtkJCNSplitField::FieldLocationNames[5][12] 
= { "DATA_OBJECT",
    "POINT_DATA",
    "CELL_DATA",
    "VERTEX_DATA",
    "EDGE_DATA" };

char vtkJCNSplitField::AttributeNames[vtkDataSetAttributes::NUM_ATTRIBUTES][10]  = { {0} };


typedef vtkJCNSplitField::Component Component;

vtkJCNSplitField::vtkJCNSplitField()
{
  this->FieldName = 0;
  this->FieldLocation = -1;
  this->AttributeType = -1;
  this->FieldType = -1;

  this->Head = 0;
  this->Tail = 0;

  //convert the attribute names to uppercase for local use
  if (vtkJCNSplitField::AttributeNames[0][0] == 0) 
    {
    for (int i = 0; i < vtkDataSetAttributes::NUM_ATTRIBUTES; i++)
      {
      int l = static_cast<int>(
        strlen(vtkDataSetAttributes::GetAttributeTypeAsString(i)));
      for (int c = 0; c < l && c < 10; c++)
        {
        vtkJCNSplitField::AttributeNames[i][c] = 
          toupper(vtkDataSetAttributes::GetAttributeTypeAsString(i)[c]);
        }
      }
    }
}

vtkJCNSplitField::~vtkJCNSplitField()
{
  delete[] this->FieldName;
  this->FieldName = 0;
  this->DeleteAllComponents();
}

void vtkJCNSplitField::SetInputField(const char* name, int fieldLoc)
{
  if (!name)
    {
    return;
    }

  if ( (fieldLoc !=  vtkJCNSplitField::DATA_OBJECT) &&
       (fieldLoc !=  vtkJCNSplitField::POINT_DATA) &&
       (fieldLoc !=  vtkJCNSplitField::CELL_DATA) &&
       (fieldLoc !=  vtkJCNSplitField::VERTEX_DATA) &&
       (fieldLoc !=  vtkJCNSplitField::EDGE_DATA) )
    {
    vtkErrorMacro("The source for the field is wrong.");
    return;
    }

  this->Modified();
  this->FieldLocation = fieldLoc;
  this->FieldType = vtkJCNSplitField::NAME;

  delete[] this->FieldName;
  this->FieldName = new char[strlen(name)+1];
  strcpy(this->FieldName, name);
}

void vtkJCNSplitField::SetInputField(int attributeType, int fieldLoc)
{
  if ( (fieldLoc !=  vtkJCNSplitField::POINT_DATA) &&
       (fieldLoc !=  vtkJCNSplitField::CELL_DATA) &&
       (fieldLoc !=  vtkJCNSplitField::VERTEX_DATA) &&
       (fieldLoc !=  vtkJCNSplitField::EDGE_DATA) )
    {
    vtkErrorMacro("The source for the field is wrong.");
    return;
    }

  this->Modified();
  this->FieldLocation = fieldLoc;
  this->FieldType = vtkJCNSplitField::ATTRIBUTE;
  this->AttributeType = attributeType;

}

void vtkJCNSplitField::SetInputField(const char* name,
                                  const char* fieldLoc)
{
  if ( !name || !fieldLoc)
    {
    return;
    }

  int numAttr = vtkDataSetAttributes::NUM_ATTRIBUTES;
  int numFieldLocs = 3;
  int i;

  // Convert strings to ints and call the appropriate SetInputField()
  int attrType=-1;
  for(i=0; i<numAttr; i++)
    {
    if (!strcmp(name, AttributeNames[i]))
      {
      attrType = i;
      break;
      }
    }

  int loc=-1;
  for(i=0; i<numFieldLocs; i++)
    {
    if (!strcmp(fieldLoc, FieldLocationNames[i]))
      {
      loc = i;
      break;
      }
    }
  if (loc == -1)
    {
    vtkErrorMacro("Location for the field is invalid.");
    return;
    }

  if (attrType == -1)
    {
    this->SetInputField(name, loc);
    }
  else
    {
    this->SetInputField(attrType, loc);
    }

}

void vtkJCNSplitField::Split(int component, const char* arrayName)
{
  if (!arrayName)
    {
    return;
    }

  this->Modified();
  Component* comp = this->FindComponent(component);
  // If component is already there, just reset the information
  if ( comp )
    {
    comp->SetName(arrayName);
    }
  // otherwise add a new one
  else
    {
    comp = new Component;
    comp->SetName(arrayName);
    comp->Index = component;
    this->AddComponent(comp);
    }
}

int vtkJCNSplitField::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);


  // get the input and output
  vtkDataObject *input = vtkDataObject::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkDataObject *output = vtkDataObject::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkDataSet *dsInput = vtkDataSet::SafeDownCast(input);
  vtkDataSet *dsOutput = vtkDataSet::SafeDownCast(output);

  vtkGraph *grInput = vtkGraph::SafeDownCast(input);
  vtkGraph *grOutput = vtkGraph::SafeDownCast(output);

  // This has to be here because it initialized all field datas.
  // output->CopyStructure( input );

  // Pass all. (data object's field data is passed by the
  // superclass after this method)
  if (dsInput)
		{
  	dsOutput->GetPointData()->PassData( dsInput->GetPointData() );
  	dsOutput->GetCellData()->PassData( dsInput->GetCellData() );
    // This has to be here because it initialized all field datas.
    dsOutput->CopyStructure( dsInput );
		}
  else if (grInput)
		{
  	grOutput->GetVertexData()->PassData( grInput->GetVertexData() );
  	grOutput->GetEdgeData()->PassData( grInput->GetEdgeData() );
    // This has to be here because it initialized all field datas.
    grOutput->CopyStructure( grInput );
		}


  
  // Pass all. (data object's field data is passed by the
  // superclass after this method)
  //output->GetPointData()->PassData( input->GetPointData() );
  //output->GetCellData()->PassData( input->GetCellData() );

  vtkDataArray* outputArray;
  vtkDataArray* inputArray = 0;
  vtkFieldData* fd = 0;
  vtkFieldData* outputFD = 0;
  Component* cur = this->GetFirst();
  Component* before;

  if (!cur) { return 1; }

  // find the input and output field data
  if ( this->FieldLocation == vtkJCNSplitField::DATA_OBJECT)
    {
    fd = input->GetFieldData();
    outputFD = output->GetFieldData();
    if (!fd || !outputFD)
      {
      vtkErrorMacro("No field data in vtkDataObject.");
      return 1;
      }
    }
  else if ( this->FieldLocation == vtkJCNSplitField::POINT_DATA )
    {
    fd = dsInput->GetPointData();
    outputFD = dsOutput->GetPointData();
    }
  else if ( this->FieldLocation == vtkJCNSplitField::CELL_DATA )
    {
    fd = dsInput->GetCellData();
    outputFD = dsOutput->GetCellData();
    }
  else if ( this->FieldLocation == vtkJCNSplitField::VERTEX_DATA )
    {
    fd = grInput->GetVertexData();
    outputFD = grOutput->GetVertexData();
    }
  else if ( this->FieldLocation == vtkJCNSplitField::EDGE_DATA )
    {
    fd = grInput->GetEdgeData();
    outputFD = grOutput->GetEdgeData();
    }

  if ( this->FieldType == vtkJCNSplitField::NAME )
    {
    inputArray = fd->GetArray(this->FieldName);
    }
  else if ( this->FieldType == vtkJCNSplitField::ATTRIBUTE )
    {
    // If we are working with attributes, we also need to have
    // access to vtkDataSetAttributes methods.
    vtkDataSetAttributes* dsa=vtkDataSetAttributes::SafeDownCast(fd);
    if (!dsa)
      {
      vtkErrorMacro("Sanity check failed, returning.");
      return 1;
      }
    inputArray = dsa->GetAttribute(this->AttributeType);
    }

  if (!inputArray)
    {
    vtkErrorMacro("Sanity check failed, returning.");
    return 1;
    }

  // iterate over all components in the linked list and 
  // generate them
  do
    {
    before = cur;
    cur = cur->Next;
    if (before->FieldName)
      {
      outputArray = this->SplitArray(inputArray, before->Index);
      if (outputArray)
        {
        outputArray->SetName(before->FieldName);
        outputFD->AddArray(outputArray);
        outputArray->UnRegister(this);
        }
      }
    } 
  while (cur);

  return 1;
}

// fast pointer copy
template <class T>
void vtkJCNSplitFieldCopyTuples(T* input, T* output, vtkIdType numTuples, 
                   int numComp, int component)
{
  for (int i=0; i<numTuples; i++)
    {
    output[i] = input[numComp*i+component];
    }
}

vtkDataArray* vtkJCNSplitField::SplitArray(vtkDataArray* da, int component)
{
  if ( (component < 0) || (component > da->GetNumberOfComponents()) )
    {
    vtkErrorMacro("Invalid component. Can not split");
    return 0;
    }

  vtkDataArray* output = da->NewInstance();
  output->SetNumberOfComponents(1);
  int numTuples = da->GetNumberOfTuples();
  output->SetNumberOfTuples(numTuples);
  if ( numTuples > 0 )
    {
    switch (output->GetDataType())
      {
      vtkTemplateMacro(
        vtkJCNSplitFieldCopyTuples((VTK_TT *)da->GetVoidPointer(0), 
                                (VTK_TT *)output->GetVoidPointer(0), 
                                numTuples,
                                da->GetNumberOfComponents(), 
                                component));
      // This is not supported by the template macro.
      // Switch to using the float interface.
      case VTK_BIT:
      {
      for(int i=0; i<numTuples; i++)
        {
        output->SetComponent(i, 0, da->GetComponent(i, component));
        }
      }
      break;
      default:
        vtkErrorMacro(<<"Sanity check failed: Unsupported data type.");
        return 0;
      }
    }
  
  return output;

}


// Linked list methods
void vtkJCNSplitField::AddComponent(Component* op)
{
  op->Next = 0;

  if (!this->Head)
    {
    this->Head = op;
    this->Tail = op;
    return;
    }
  this->Tail->Next = op;
  this->Tail = op;
}

Component* vtkJCNSplitField::FindComponent(int index)
{
  Component* cur = this->GetFirst();
  if (!cur) { return 0; }

  if (cur->Index == index) { return cur; }
  while (cur->Next)
    {
    if (cur->Next->Index == index)
      {
      return cur->Next;
      }
    cur = cur->Next;
    }
  return 0;
}

void vtkJCNSplitField::DeleteAllComponents()
{
  Component* cur = this->GetFirst();
  if (!cur) {return;}
  Component* before;
  do 
    {
    before = cur;
    cur = cur->Next;
    delete before;
    } 
  while (cur);
  this->Head = 0;
  this->Tail = 0;
}

void vtkJCNSplitField::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "Field name: ";
  if (this->FieldName)
    {
    os << this->FieldName << endl;
    }
  else
    {
    os << "(none)" << endl;
    }
  os << indent << "Field type: " << this->FieldType << endl;
  os << indent << "Attribute type: " << this->AttributeType << endl;
  os << indent << "Field location: " << this->FieldLocation << endl;
  os << indent << "Linked list head: " << this->Head << endl;
  os << indent << "Linked list tail: " << this->Tail << endl;
  os << indent << "Components: " << endl;
  this->PrintAllComponents(os, indent.GetNextIndent());
}

void vtkJCNSplitField::PrintComponent(Component* op, ostream& os,
                                   vtkIndent indent)
{
  os << indent << "Field name: " << op->FieldName << endl;
  os << indent << "Component index: " << op->Index << endl;
}

void vtkJCNSplitField::PrintAllComponents(ostream& os, vtkIndent indent)
{
  Component* cur = this->GetFirst();
  if (!cur) { return; }
  Component* before;
  do
    {
    before = cur;
    cur = cur->Next;
    os << endl;
    this->PrintComponent(before, os, indent);
    } 
  while (cur);
}
