/*=========================================================================

 *	File: vtkRAWReader.cxx
 *	Graph visualization library for VTK
 * 
 *	This software is distributed WITHOUT ANY WARRANTY; 
 *	without even the implied warranty of MERCHANTABILITY 
 *	or FITNESS FOR A PARTICULAR PURPOSE.  
 * 
 *	See the file copyright.txt for details.  

=========================================================================*/

#include "vtkErrorCode.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkIdList.h"
#include "vtkInformationVector.h"
#include "vtkRAWReader.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredPoints.h"
#include "vtkSmartPointer.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkUnsignedCharArray.h"
#include <vtksys/ios/sstream>
#include <sys/stat.h>
#include <vector>

typedef std::vector<vtkSmartPointer<vtkIdList> > *Polytope;

vtkStandardNewMacro(vtkRAWReader);

vtkRAWReader::vtkRAWReader()
{
  this->SetNumberOfInputPorts(0);
  this->FileName = NULL;
  this->DimX = 0;
  this->DimY = 0;
  this->DimZ = 0;
}

vtkRAWReader::~vtkRAWReader()
{
  if (this->FileName)
    {
    delete [] this->FileName;
    }
}

void vtkRAWReader::SetRAWFile(char *fileName, int dimX, int dimY, int dimZ)
{
  this->SetFileName(fileName);
  this->SetDimX(dimX);
  this->SetDimY(dimY);
  this->SetDimZ(dimZ);
}

//----------------------------------------------------------------------------
// Default method performs Update to get information.  Not all the old
// structured points sources compute information

int vtkRAWReader::RequestInformation (
  vtkInformation * vtkNotUsed(request),
  vtkInformationVector ** vtkNotUsed( inputVector ),
  vtkInformationVector *outputVector)
{
  int extent[6];
  
  extent[0] = 0; extent[1] = this->DimX - 1;
  extent[2] = 0; extent[3] = this->DimY - 1;
  extent[4] = 0; extent[5] = this->DimZ - 1;
  
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  outInfo->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);

  this->GetOutput()->SetNumberOfScalarComponents(1,outInfo);
  this->GetOutput()->SetScalarType(VTK_UNSIGNED_CHAR,outInfo);

  return 1;
}


//-----------------------------------------------------------------------------
int vtkRAWReader::RequestData(
  vtkInformation *,
  vtkInformationVector **,
  vtkInformationVector *outputVector)
{
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkImageData *output = vtkImageData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  output->SetExtent(
    outInfo->Get(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT()));
  output->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

  std::vector<Polytope> pointers 
      = std::vector<Polytope>();
      
  std::vector<vtkIdList*> dumb 
      = std::vector<vtkIdList*>();
      
  Polytope frag;
        
  vtkIdList *dlist;      
  vtkSmartPointer<vtkIdList> list;

  if (!this->DimZ || !this->DimY || !this->DimZ)
    {
    vtkErrorMacro("Dataset dimension not set.");
    return 0;
    }
  if (!this->FileName || (strlen(this->FileName) == 0))
    {
    vtkErrorMacro(<< "A FileName must be specified.");
    this->SetErrorCode(vtkErrorCode::NoFileNameError);
    return 0;
    }
    
  struct stat fs;
  istream *is = NULL;
  
  if (stat(this->FileName, &fs) != 0)
    {
    vtkErrorMacro(<< "Unable to open file: "<< this->FileName);
    this->SetErrorCode( vtkErrorCode::CannotOpenFileError );
    return 0;
    }

  is = new ifstream(this->FileName, ios::in);
  if (is->fail())
    {
    vtkErrorMacro(<< "Unable to open file: "<< this->FileName);
    delete is;
    is = NULL;
    this->SetErrorCode( vtkErrorCode::CannotOpenFileError );
    return 0;
    }

  int numValues = this->DimX * this->DimY * this->DimZ;
  
  vtkSmartPointer<vtkUnsignedCharArray> array 
      = vtkSmartPointer<vtkUnsignedCharArray>::New();
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(numValues);
  
  unsigned char *ptr = ((vtkUnsignedCharArray *)array)->WritePointer(0,numValues);
 
  is->read((char *)ptr, sizeof(unsigned char)*numValues);
  
  if (is->eof())
    {
    vtkGenericWarningMacro(<<"Error reading binary data!");
    delete is;
    return 0;
    }
  vtkDebugMacro (<< "Finished reading RAW file");
  delete is;
  
  output->GetPointData()->SetScalars(array);
  output->GetPointData()->GetScalars()->SetName("RAW");
  
  for (int i = 0; i < this->DimX; i++)
    {
    frag = new std::vector<vtkSmartPointer<vtkIdList> >();
    pointers.push_back(frag);
    for (int j = 0; j < this->DimY; j++)
      {
      list = vtkSmartPointer<vtkIdList>::New();
      list->SetNumberOfIds(10);
      frag->push_back(list);
      pointers[0]->push_back(list);
      }
    }    

  for (std::vector<Polytope>::iterator it = pointers.begin();
      it != pointers.end();
      it++)
    {
    frag = *it;
    for (std::vector<vtkSmartPointer<vtkIdList> >::iterator fc = frag->begin();
        fc != frag->end();
        fc++)
      {
      (*fc)->Delete();
      }
    }    
  return 1;
}


void vtkRAWReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "File Name: "
     << (this->FileName ? this->FileName : "(none)") << endl;
  os << indent << "Extent: "
     << this->DimX << " X " << this->DimY << " X " << this->DimZ
     << endl; 
}

