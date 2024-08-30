#Reference:
#[1] Agarwal, T., Chattopadhyay, A., Natarajan, V. (2021).
#Topological Feature Search in Time-Varying Multifield Data.
#In: Topological Methods in Data Analysis and Visualization VI. Mathematics and Visualization. Springer, Cham.

from vtkMETAPython import vtkComputeUnstructuredGrid, vtkJointContourNetForShapes, vtkMultiDimensionalReebGraph
import vtk
from math import *
from sys import exit
from sys import argv
import re


#Parameters
fields=['NormalizedGeodesicDistance','NormalizedD2']
shapeFileName1 = 'Shape1'
shapeFileName2 = 'Shape2'
numberOfSlabs = 32
w = 1.0
q = 1.0



#Read the function values at the vertices of a shape
def ReadFieldValues(fileName):
	valueList=[]
	for i in range(0,len(fields)):
		fieldFileName =  '/app/data/Shapes/'+ fileName+'_'+fields[i]+'.txt'
		fieldFile = open(fieldFileName)
		lines = fieldFile.readlines()
		values=[]
		for line in lines:
			value = line.rstrip('\n')
			values.append(round(float(value)*pow(10,4),3))
		valueList.append(values)
		fieldFile.close()
	return valueList

#Construct the simiplical complex
def ConstructSimplicialComplex(fileName, fieldValues):
	file1 = open('/app/data/Shapes/' + fileName)
	file1.readline()
	entries = file1.readline().rstrip('\n').split(' ')
	points = int(entries[0])
	numberOftriangles = int(entries[1])
	pointList = vtk.vtkFloatArray()
	for i in range(0,points):
		point  = file1.readline().rstrip('\n').split(' ')
		point[0] = float(point[0])
		point[1] = float(point[1])
		point[2] = float(point[2])
		pointList.InsertNextValue(point[0])
		pointList.InsertNextValue(point[1])
		pointList.InsertNextValue(point[2])
	triangleList = vtk.vtkIdTypeArray()
	for i in range(0,numberOftriangles):
		vertexIDs = file1.readline().rstrip('\n').split(' ')
		if(vertexIDs[0]!='3'):
			print('Error in File')
			exit()
		vertex1 = int(float(vertexIDs[1]))
		vertex2 = int(float(vertexIDs[2]))
		vertex3 = int(float(vertexIDs[3]))
		triangleList.InsertNextValue(vertex1)
		triangleList.InsertNextValue(vertex2)
		triangleList.InsertNextValue(vertex3)
	file1.close()							
	pointData = vtk.vtkPointData()
	for i in range(0,len(fields)):
		f = vtk.vtkFloatArray()
		f.SetName("f"+str(i+1))
		for p in range(0, points):
			f.InsertNextValue(fieldValues[i][p])
		pointData.AddArray(f)
	simp = vtkComputeUnstructuredGrid()
	simp.SetInputData(vtk.vtkImageData())
	simp.Update()
	simp.ComputeUnstructuredGridFromTriangles(pointList, pointData,triangleList)
	return simp

#Computing the number of quantized fiber-components in each of the range-value bins
def ComputeRangeBinCount(jcnGraph):
	range_bin_count={}
	for i in range(jcnGraph.GetNumberOfVertices()):
		valueList=[]
		for f in range(len(fields)):
			value = jcnGraph.GetVertexData().GetArray("f"+str(f+1)).GetComponent(i,0)
			valueList.append(value)
		if(tuple(valueList) not in range_bin_count):
				range_bin_count[tuple(valueList)] = 1.0
		else:
				range_bin_count[tuple(valueList)] = range_bin_count[tuple(valueList)] +1.0
	return range_bin_count


if __name__ == "__main__":

	shape1FieldValues = ReadFieldValues(shapeFileName1)
	shape2FieldValues = ReadFieldValues(shapeFileName2)

	#Constructing simplicial complexes
	simplicialComplex1 = ConstructSimplicialComplex(shapeFileName1 + '.off', shape1FieldValues)
	simplicialComplex2 = ConstructSimplicialComplex(shapeFileName2 + '.off', shape2FieldValues)

	#Computing JCNs
	jcn1 = vtkJointContourNetForShapes()
	jcn1.SetInputConnection(simplicialComplex1.GetOutputPort())

	jcn2 = vtkJointContourNetForShapes()
	jcn2.SetInputConnection(simplicialComplex2.GetOutputPort())

	for f in range(len(fields)):
		minValue = floor(min(min(shape1FieldValues[f]), min(shape2FieldValues[f])))-1
		maxValue = ceil(max(max(shape1FieldValues[f]), max(shape2FieldValues[f])))+1
		slabWidth = float(maxValue-minValue)/numberOfSlabs
		fieldName = 'f'+str(f+1)
		jcn1.AddField(fieldName, slabWidth, minValue)
		jcn2.AddField(fieldName, slabWidth, minValue)

	jcn1.Update()
	jcn2.Update()

	jcnGraph1 = jcn1.GetOutput(0)
	jcnGraph2 = jcn2.GetOutput(0)


	#Computing MDRGs
	mdrg1 =  vtkMultiDimensionalReebGraph() 
	mdrg2 =  vtkMultiDimensionalReebGraph() 
	mdrg1.SetInputData(jcnGraph1)
	mdrg2.SetInputData(jcnGraph2)
	for f in range(len(fields)):
		fieldName = 'f'+str(f+1)
		mdrg1.AddField(fieldName)
		mdrg2.AddField(fieldName)
	mdrg1.Update()
	mdrg2.Update()

	#Computing Histograms corresponding to the fiber component distributions
	range_bin_count_jcn_1=ComputeRangeBinCount(jcnGraph1)
	range_bin_count_jcn_2=ComputeRangeBinCount(jcnGraph2)

	for value in range_bin_count_jcn_1:
		if(value not in range_bin_count_jcn_2):
			range_bin_count_jcn_2[value]=0

	for value in range_bin_count_jcn_2:
		if(value not in range_bin_count_jcn_1):
			range_bin_count_jcn_1[value]=0

	#Computing range values corresponding to jacobi nodes
	jacobi_range_values=set()
	for i in range(jcnGraph1.GetNumberOfVertices()):
		if(mdrg1.GetOutput().GetVertexData().GetArray("MDRG").GetValue(i)==1):
			valueList=[]
			for f in range(len(fields)):
				value = jcnGraph1.GetVertexData().GetArray("f"+str(f+1)).GetComponent(i,0)
				valueList.append(value)
			jacobi_range_values.add(tuple(valueList))

	for i in range(jcnGraph2.GetNumberOfVertices()):
		if(mdrg2.GetOutput().GetVertexData().GetArray("MDRG").GetValue(i)==1):
			valueList=[]
			for f in range(len(fields)):
				value = jcnGraph2.GetVertexData().GetArray("f"+str(f+1)).GetComponent(i,0)
				valueList.append(value)
			jacobi_range_values.add(tuple(valueList))

	#Computing the Wasserstein distance (see Equation (14) in [1])
	sum = 0
	for value in range_bin_count_jcn_1:
		if(value in jacobi_range_values):
			sum = sum + w*(pow(abs( (range_bin_count_jcn_1[value]/jcnGraph1.GetNumberOfVertices()) - (range_bin_count_jcn_2[value]/jcnGraph2.GetNumberOfVertices()) ), q))
		else:
			sum = sum + pow(abs( (range_bin_count_jcn_1[value]/jcnGraph1.GetNumberOfVertices()) - (range_bin_count_jcn_2[value]/jcnGraph2.GetNumberOfVertices()) ), q)

	distance = pow(sum,1.0/q)
	# print(str('Distance between Shape1 and Shape2: ') + str(distance) + '\n')
	print(distance)
