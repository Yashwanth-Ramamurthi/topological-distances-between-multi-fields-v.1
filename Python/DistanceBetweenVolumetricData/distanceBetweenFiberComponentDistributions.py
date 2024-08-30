#Reference:
#[1] Agarwal, T., Chattopadhyay, A., Natarajan, V. (2021).
#Topological Feature Search in Time-Varying Multifield Data.
#In: Topological Methods in Data Analysis and Visualization VI. Mathematics and Visualization. Springer, Cham.


from vtkMETAPython import vtkSimplicate, vtkJointContourNet, vtkRAWReader, vtkMultiDimensionalReebGraph
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
numberOfSlabs = 16

w = 1.0
q = 1.0


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

	#Adding the parameters for constructing the JCNs

	for f in range(len(fields)):
		minValue = floor(min(min(data1_FieldValues[f]), min(data2_FieldValues[f])))-1
		maxValue = ceil(max(max(data1_FieldValues[f]), max(data2_FieldValues[f])))+1
		slabWidth=float(maxValue-minValue)/numberOfSlabs	
		jcn1.AddField(fieldNames[f], slabWidth, minValue)
		jcn2.AddField(fieldNames[f], slabWidth, minValue)

	#Computing the JCNs
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

	# print(str('Distance between the two bivariate data: ') + str(distance) + '\n')

	print(distance)
