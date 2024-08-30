#Reference:
#[1] Yashwanth Ramamurthi, Tripti Agarwal and Amit Chattopadhyay. "A Topological Similarity Measure Between Multi-Resolution Reeb Spaces," in IEEE Transactions on Visualization and Computer Graphics, vol. 28, no. 12, pp. 4360-4374, 1 Dec. 2022


from vtkMETAPython import vtkSimplicate, vtkJointContourNet, vtkSimilarityBetweenMultiResolutionReebSpaces, vtkRAWReader
import vtk
from math import *
from sys import exit
from sys import argv
import re


#Parameters
data_dimension = [19,19,19]
fieldNames = ['f1', 'f2']
timesteps = ['260','270']
fields = ['p','n']

weightVolume = 0.25
weightRangeMeasure = 0.25
weightBettiZero = 0.25
weightDegree = 0.25

numberOfResolutions = 5



#Read the function values at the vertices of a shape
def ReadFieldValues(timestep):
	pointData = vtk.vtkImageData()
	pointData.SetDimensions(data_dimension[0], data_dimension[1], data_dimension[2])

	for f in range(len(fieldNames)):
		raw = vtkRAWReader()
		raw.SetDimX(data_dimension[0])
		raw.SetDimY(data_dimension[1])
		raw.SetDimZ(data_dimension[2])
		raw.SetFileName('/app/data/VolumetricData/fm258skm-Sc-00' + timestep+'_rho_'+fields[f]+'.raw')
		raw.Update()

		field_values = raw.GetOutput().GetPointData().GetArray("RAW")
		field_values.SetName(fieldNames[f])
		pointData.GetPointData().AddArray(field_values)
	return pointData

#Construct the simiplical complex
def ConstructSimplicialComplex(pointData):
	simplicialComplex = vtkSimplicate()
	simplicialComplex.SetInputData(pointData)
	simplicialComplex.SetScheme(0)
	simplicialComplex.Update()
	return simplicialComplex

if __name__ == "__main__":
	#Reading the field values
	pointData1 = ReadFieldValues(timesteps[0])
	pointData2 = ReadFieldValues(timesteps[1])

	n = pointData2.GetPointData().GetArray("f1").GetNumberOfTuples()

	#Storing the field values
	data1_FieldValues=[]
	data2_FieldValues=[]
	for f in range(len(fieldNames)):
		values=[]
		for i in range(n):
			currentValue = pointData1.GetPointData().GetArray(fieldNames[f]).GetComponent(i,0)
			values.append(currentValue)
		data1_FieldValues.append(values)

	for f in range(len(fieldNames)):
		values=[]
		for i in range(n):
			currentValue = pointData2.GetPointData().GetArray(fieldNames[f]).GetComponent(i,0)
			values.append(currentValue)
		data2_FieldValues.append(values)

	#Constructing simplicial complexes
	simplicialComplex1 = ConstructSimplicialComplex(pointData1)
	simplicialComplex2 = ConstructSimplicialComplex(pointData2)

	#Computing JCNs
	jcn1 = vtkJointContourNet()
	jcn1.SetInputConnection(simplicialComplex1.GetOutputPort())

	jcn2 = vtkJointContourNet()
	jcn2.SetInputConnection(simplicialComplex2.GetOutputPort())

	#Adding the parameters for constructing the JCN, computing the distance between MRSs
	similarityBetweenMRSs = vtkSimilarityBetweenMultiResolutionReebSpaces()
	similarityBetweenMRSs.SetAttributes(weightVolume, weightRangeMeasure, weightBettiZero, weightDegree)

	for f in range(len(fields)):
		minValue = floor(min(min(data1_FieldValues[f]), min(data2_FieldValues[f])))-1
		maxValue = ceil(max(max(data1_FieldValues[f]), max(data2_FieldValues[f])))+1
		slabWidth=float(maxValue-minValue)/pow(2.0,numberOfResolutions-1)
		jcn1.AddField(fieldNames[f], slabWidth, minValue)
		jcn2.AddField(fieldNames[f], slabWidth, minValue)
		similarityBetweenMRSs.AddParameters(fieldNames[f], minValue, maxValue, slabWidth)

	#Computing the JCNs
	jcn1.Update()
	jcn2.Update()

	jcnGraph1 = jcn1.GetOutput(0)
	jcnGraph2 = jcn2.GetOutput(0)

	#Computing the similarity between MRSs
	similarity = similarityBetweenMRSs.ComputeSimilarity(jcnGraph1, jcnGraph2, numberOfResolutions)


	#Distance between MRSs (see Equation (7) in [1])
	distance = 1 - similarity
	# print(str('Distance between the two bivariate data: ') + str(distance) + '\n')
	print(distance)
