/*=========================================================================
This source has no copyright.  It is intended to be copied by users
wishing to create their own VTK classes locally.
=========================================================================*/

#include "vtkObjectFactory.h"
#include "vtkPyramidTree.h"

vtkStandardNewMacro(vtkPyramidTree);

//----------------------------------------------------------------------------

vtkPyramidTree::vtkPyramidTree()
{
  this->ColNumber = 0;
  this->TTime = 0;
  this->Interval = 1;
}

//----------------------------------------------------------------------------

void vtkPyramidTree::Initialize(int col, const double *bounds)
{
  this->RawData.clear();
  this->SumData.clear();
  this->Map.clear();
  if (col != this->ColNumber)
    {
    this->ColNumber = col;
    this->MaxDim.resize(this->ColNumber);
    this->MinDim.resize(this->ColNumber);
    this->Median.resize(this->ColNumber);
    double myInterval = static_cast<double>(MAX_BUCKET_NUMBER/(this->ColNumber*2));
    this->Interval = floor(myInterval);
    double m = static_cast<double>(1.0 / this->ColNumber);
    int div = ceil(pow(MAX_BUCKET_NUMBER,m));
    for (int i = 0; i < this->ColNumber; i++)
      {
      this->Median[i] = div;
      }
    }
    
  for (int i = 0; i < this->ColNumber; i++)
     {
      this->MinDim[i] = bounds[2*i];
      this->MaxDim[i] = bounds[2*i+1];
      if (MaxDim[i] <= MinDim[i])
        {
        MaxDim[i] = MinDim[i] + 1;
        }
     }
}

//----------------------------------------------------------------------------
vtkPyramidTree::~vtkPyramidTree()
{
}

//----------------------------------------------------------------------------
void vtkPyramidTree::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//------------------
int vtkPyramidTree::GetNumberElement()
{
  return this->RawData.size();
}

//-----------------
int vtkPyramidTree::GetNumberDimension()
{
  return this->ColNumber;
}

//----------------
void vtkPyramidTree::SetDimension(int dim)
{
  this->ColNumber = dim;
}

//------------------
void vtkPyramidTree::SetMaxBounds(const double *bounds)
{
  for (int i = 0; i < this->ColNumber; i++)
     {
      this->MaxDim[i] = bounds[i];
     }
}

//------------------
void vtkPyramidTree::SetMinBounds(const double *bounds)
{
  for (int i = 0; i < this->ColNumber; i++)
    {
    this->MinDim[i] = bounds[i];
    }
}

//------------------
int vtkPyramidTree::InsertUniquePoint(const double *point,vtkIdType &pid)
{
  std::vector<double> myPoint(this->ColNumber);
    
  for (int i = 0; i < this->ColNumber; i++)
    {
    myPoint[i] = point[i];
    }

  vtkIdType id;
  this->GetPointIndex(id, point);
  if (id < 0)
    {
    this->InsertToStructure(myPoint,id);
    pid = this->RawData.size() -1 ;
    return 1;
    }
  else
    {
    pid = id;
    return 0;
    }
}

//-------------------------

void vtkPyramidTree::InsertToStructure(std::vector<double> point, vtkIdType duplicate)
{
  // Transform the original data to the encoded data sets.
  double sum = ComputePTStruture(point);
    
  this->RawData.push_back(point);
  this->SumData.push_back(sum);
    
  int searchKey = this->HashValue(point);
  int currentIndex = this->RawData.size() - 1;
    
  if (duplicate == -1)
    {
    // If there is no inserted data, have not search key
    std::vector<int> myVector;
    myVector.push_back(currentIndex);
    this->Map[searchKey] = myVector;
    }
  else if (duplicate == -2)
    {
    //If there is no inserted data, have search key
    std::vector<int> value = this->Map.find(searchKey)->second;
    value.push_back(currentIndex);
    this->Map[searchKey] = value;
    }
}


//------------------
void vtkPyramidTree::GetPointIndex(vtkIdType &pid, const double *point)
{
  int finalIndex;
  std::vector<double> searchValue(this->ColNumber);
  for (int i = 0; i < this->ColNumber; i++)
    {
    searchValue[i]= point[i];
    }
    
  double mySum = ComputePTStruture(searchValue);  
  int searchKey = this->HashValue(searchValue);
  if (this->Map.find(searchKey) != this->Map.end())
    {
    std::vector<int> output = this->Map.find(searchKey)->second;
    bool contain = false;
    
    for (int i = 0; i < output.size(); i++)
      {
      int myIndex = output[i];
      double sum = this->SumData[myIndex];
      if (!compare(sum,mySum))
        {
        std::vector<double> myPoint = this->RawData[myIndex];
        if (!compareVector(searchValue,myPoint))
          {
          // if there is such data
          finalIndex = myIndex;
          contain = true;
          break;
          }
        }
      }
    // if there is no inserted data, but there is search key
    if (!contain) 
      {
      finalIndex = -2;
      }
    }
  else
    {
    // if there is no inserted data
    finalIndex = -1;
    }
  pid = finalIndex;
}


//----------------------------------------------------------------------------
void vtkPyramidTree::GetPoint(vtkIdType pid, double *point)
{
  int size = this->RawData.size();
  if (pid < 0 || pid > size)
    {
    vtkErrorMacro("ID out of bounds");
    exit(1);
    }
  std::vector<double> my = this->RawData[pid];
  for (int i = 0; i < this->ColNumber; i++)
     {
      point[i] = my[i];
     }
}

//----------------------------------------------------------------------------

int vtkPyramidTree::GetNumberOfPoints()
{
  return this->RawData.size();
}

//----------------------------------------------------------------------------
void vtkPyramidTree::InsertNextPoint(const double *point)
{
  std::vector<double> temp;
  for (int i = 0; i < this->ColNumber; i++)
    {
    temp.push_back(point[i]);
    }
  this->RawData.push_back(temp);
}

//----------------------------------------------------------------------------
double vtkPyramidTree::GetTime()
{
  return TTime;
}

//----------------------------------------------------------------------------
int vtkPyramidTree::ComputePTStruture(std::vector<double> myPoint)
{
  // Transform the original data to the encoded data sets.
  double tempValue[this->ColNumber];
  for (int j = 0; j < this->ColNumber; j++)
    {
    double ele = myPoint[j];
    ele = (ele - this->MinDim[j]) / (this->MaxDim[j] - this->MinDim[j]);
    tempValue[j] = ele;
    }
    
  // Calculate the interpolted 1-D point
  int dMax = 0;
  double diff = 0.5 -tempValue[0];
  double height = abs(diff);
  int index;
  for (int j = 1; j < this->ColNumber; j++)
    {
    double diff = 0.5-tempValue[j];
    if (height < abs(diff))
      {
      dMax = j;
      height = abs(diff);
      }
    }
    
  if (tempValue[dMax] < 0.5)
    {
    index = dMax;
    }
  else 
    {
    index = dMax + this->ColNumber;
    }
    
  index = index * this->Interval;
  return index;  
}

//----------------------------------------------------------------------------
int vtkPyramidTree::HashValue(std::vector<double> point)
{
  int searchKey = 0 ;
  int value[this->ColNumber];
  for (int i = 0; i < this->ColNumber; i++ )
    {
    value[i] = static_cast<int>(
        static_cast<double>((point[i] - this->MinDim[i]) /
                            (this->MaxDim[i] - this->MinDim[i]))
            * this->Median[i]);
    if (value[i] >= this->Median[i])
      {
      value[i] = this->Median[i] - 1;
      }
    }
    
  for (int i = 0; i < this->ColNumber; i++)
    {
    int temp = value[i];
    for (int j = 0; j < i; j++)
      {
      temp = temp * this->Median[j];
      }
    searchKey = searchKey + temp;
    }
  return searchKey;
}

int vtkPyramidTree::GetBucketCount()
{
  return this->Map.bucket_count();  
}

int vtkPyramidTree::GetBucketSize()
{
  int sum = 0;
  for (int i = 0; i < this->GetBucketCount(); i++)
     {
      sum++;
     }
  return sum;
}
