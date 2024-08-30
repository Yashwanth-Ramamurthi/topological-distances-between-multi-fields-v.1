#Reference:
#[1] Yashwanth Ramamurthi and Amit Chattopadhyay. 2023. Topological Shape Matching using Multi-Dimensional Reeb Graphs. In Proceedings of the Thirteenth Indian Conference on Computer Vision, Graphics and Image Processing (ICVGIP '22).
#####Association for Computing Machinery, New York, NY, USA, Article 5, 1-10.


from vtkMETAPython import vtkSimplicate, vtkJointContourNet, vtkDistanceBetweenMultiDimensionalReebGraphs, vtkRAWReader
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

w_0 = 1.0/3
w_1 = 1.0/3
w_2 = 1 - (w_0 + w_1)

numberOfSlabs = 16

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

	#Adding the parameters for constructing the JCN, computing the distance between MDRGs
	#Computation of the three distances in Equation (11) in [1]
	distanceBetweenMDRGs_sublevelset =  vtkDistanceBetweenMultiDimensionalReebGraphs() 
	distanceBetweenMDRGs_superlevelset =  vtkDistanceBetweenMultiDimensionalReebGraphs() 
	distanceBetweenMDRGs_extended =  vtkDistanceBetweenMultiDimensionalReebGraphs()

	#Adding the parameters for constructing the JCN, computing the distance between MRSs

	for f in range(len(fields)):
		minValue = floor(min(min(data1_FieldValues[f]), min(data2_FieldValues[f])))-1
		maxValue = ceil(max(max(data1_FieldValues[f]), max(data2_FieldValues[f])))+1
		slabWidth=float(maxValue-minValue)/numberOfSlabs	
		jcn1.AddField(fieldNames[f], slabWidth, minValue)
		jcn2.AddField(fieldNames[f], slabWidth, minValue)
		distanceBetweenMDRGs_sublevelset.AddField(fieldNames[f])
		distanceBetweenMDRGs_superlevelset.AddField(fieldNames[f])
		distanceBetweenMDRGs_extended.AddField(fieldNames[f])

	#Computing the JCNs
	jcn1.Update()
	jcn2.Update()

	jcnGraph1 = jcn1.GetOutput(0)
	jcnGraph2 = jcn2.GetOutput(0)

	#Computation of the three distances in Equation (11) in [1]
	distance1 = distanceBetweenMDRGs_sublevelset.ComputeDistance(jcnGraph1, jcnGraph2, 0)
	distance2 = distanceBetweenMDRGs_superlevelset.ComputeDistance(jcnGraph1, jcnGraph2, 1)
	distance3 = distanceBetweenMDRGs_extended.ComputeDistance(jcnGraph1, jcnGraph2, 3)

	totalDistance = (w_0*distance1) + (w_1*distance2) + (w_2*distance3)

	# print(str('Distance between the two bivariate data: ') + str(totalDistance) + '\n')
	print(totalDistance)
