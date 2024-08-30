#Reference:
#[1] Yashwanth Ramamurthi and Amit Chattopadhyay, "A Topological Distance Between Multi-Fields Based on Multi-Dimensional Persistence Diagrams," in IEEE Transactions on Visualization and Computer
#Graphics, vol. 30, no. 9, pp. 5939-5952, Sept. 2024


from vtkMETAPython import vtkComputeUnstructuredGrid, vtkJointContourNetForShapes, vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs
import vtk
from math import *
from sys import exit
from sys import argv
import re



#Parameters
fields=['NormalizedGeodesicDistance','NormalizedD2']
shapeFileName1 = 'Shape1'
shapeFileName2 = 'Shape2'
q = 1.0
numberOfSlabs = 32

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


	wassersteinDistance =  vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs()

	#Adding the parameters for constructing the JCNs, computing the Wasserstein distance between MDRGs
	for f in range(len(fields)):
		minValue = floor(min(min(shape1FieldValues[f]), min(shape2FieldValues[f])))-1
		maxValue = ceil(max(max(shape1FieldValues[f]), max(shape2FieldValues[f])))+1
		slabWidth = float(maxValue-minValue)/numberOfSlabs
		fieldName = 'f'+str(f+1)
		jcn1.AddField(fieldName, slabWidth, minValue)
		jcn2.AddField(fieldName, slabWidth, minValue)
		wassersteinDistance.AddField(fieldName)

	wassersteinDistance.SetParameterQ(q)

	#Computing the JCNs
	jcn1.Update()
	jcn2.Update()

	jcnGraph1 = jcn1.GetOutput(0)
	jcnGraph2 = jcn2.GetOutput(0)

	#Computing the Wasserstein distance between MDRGs
	distance = wassersteinDistance.ComputeDistance(jcnGraph1, jcnGraph2)
	# print('Distance between Shape1 and Shape2: '+ str(distance) + '\n')
	print(distance)