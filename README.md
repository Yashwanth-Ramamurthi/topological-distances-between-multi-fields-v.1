File README.md

##Background - Topological distances between multi-fields

This archive contains source code and demonstration scripts for 
computing the following topological distances between multi-fields based on 
the corresponding quantized Reeb spaces (or Joint Contour Nets):
	(i) Distance between fiber-component distributions [1]
	(ii) Distance between multi-resolution Reeb spaces [2]
	(iii) Distance between multi-dimensional Reeb graphs [3]
	(iv) Distance between multi-dimensional persistence diagrams [4]

##Prerequistes
#Installation of Docker
Docker needs to be installed for running the codes in this software
Please see the following webpage for docker installation: https://docs.docker.com/engine/install/

#Memory Requirement
The installation of this software needs additional softwares/libraries, which together require 3.5 GB of disk memory.
For more details, please refer to the 'Dockerfile'.

##Installation
1. Start the Docker Engine
2. Open the Command Prompt / Terminal
3. Go to the folder 'Docker-Topological-distances-between-multi-fields' (the folder containing this file)
4. Execute the command 'docker build -t topological-distances .' (this will build the docker image)


##Running the codes
1. Start the Docker Engine (if not started already)
2. Open the Command Prompt / Terminal
3. Computing distance between shape data: Execute the command 'docker run -it --rm topological-distances shape SCRIPT_NAME'
4. Computing distance between volumetric data: Execute the command 'docker run -it --rm topological-distances volumetric SCRIPT_NAME' 
5. SCRIPT_NAME can be one of the following:
	(i) distanceBetweenFiberComponentDistributions [1]
	(ii) distanceBetweenMRSs [2]
	(iii) distanceBetweenMDRGs [3]
	(iv) distanceBetweenMDPDs [4]

##Test Data Description

#Shape Data
1. The 'TestData/Shapes' consists of two 3D shapes (in OFF format).
2. The functions normalized geodesic distance and normalized Euclidean distance (D2) are computed on each of the shapes (stored in text files)
3. The python scripts in 'Python/DistanceBetweenShapes' compute distances between the two shapes based on the bivariate field consisting of the functions mentioned above.
4. Please see Section 5.2 in [4] for more details on using these functions as shape descriptors

#Volumetric Data
1. The 'TestData/VolumetricData' folder consists of data corresponding to two timesteps of the time-varying data of the Fermium-258 atom
2. The proton, neutron, and nucleon (total) densities of the Fermium-258 atom for each of these timesteps are stored in separate files.
3. The 'Python/DistanceBetweenVolumetricData' folder contains the python scripts for computing distances between the data at the two timesteps, using the bivariate field consisting of the proton and neutron densities
4. Please see Section 6.3 in [1] for more details on the data

##Differences between the computations of the simplicial complexes or JCNs for shape and volumetric data:
	(i) The construction of similicial complexes for (2D or 3D) volumetric data and shape data, are performed by different classes (vtkSimplicate and vtkComputeUnstructuredGrid).
	(ii) The construction of joint contour nets (JCNs) for (2D or 3D) volumetric data and shape data, are performed by different classes (vtkJointContourNet and vtkJointContourNetForShapes).
Please see the python scripts in the folders 'Python/DistanceBetweenShapes' and 'Python/DistanceBetweenVolumetricData' for more details.

##Cuztomizing the code for other data:
1. Create a docker volume by executing the command 'docker volume create distances-data'
2. Copy the TestData into the volume by running this command 'cp -R path/to/your/dataset/* /var/lib/docker/volumes/distances-data/_data ' eg:  'cp -R TestData/* /var/lib/docker/volumes/distances-data/_data'. (only for linux. Needs to be modified for other operating systems)

3. Change the python scripts, the DockerFile, and the shell script 'compute-distance.sh' accordingly  
4. Build the docker: Execute the command 'docker build -t topological-distances .'
5. Run the command 'docker run -it --rm -v distances-data:/app/data topological-distances shape SCRIPT_NAME'. 

##Thanks

Prof. Hamish Carr and Prof. David Duke, University of Leeds, United Kingdom,
    for providing the software required for computing the Joint
    Contour Net.  The current software has been built on top of the software for computing the Joint Contour Net, 
    developed as a part of the project "Multifield Extension of
    Topological Analysis (META)".

##Acknowledgements
This work was funded by the Science and Engineering Research Board (SERB), India,
Grant Nr, SERB/CRG/2018/000702

Thanks to International Institute of Information Technology (IIITB),
Bangalore for funding this work. 

##References
If the distances are helpful for your research, please cite the corresponding papers.

[1] Tripti Agarwal, Amit Chattopadhyay and Natarajan (2021). Topological Feature Search in Time-Varying Multifield Data. In: Topological Methods in Data Analysis and Visualization VI. Mathematics and Visualization. Springer, Cham.
[2] Yashwanth Ramamurthi, Tripti Agarwal and Amit Chattopadhyay. "A Topological Similarity Measure Between Multi-Resolution Reeb Spaces," in IEEE Transactions on Visualization and Computer Graphics, vol. 28, no. 12, pp. 4360-4374, 1 Dec. 2022
[3] Yashwanth Ramamurthi and Amit Chattopadhyay. 2023. Topological Shape Matching using Multi-Dimensional Reeb Graphs. In Proceedings of the Thirteenth Indian Conference on Computer Vision, Graphics and Image Processing (ICVGIP '22). Association for Computing Machinery, New York, NY, USA, Article 5, 1-10.
[4] Yashwanth Ramamurthi and Amit Chattopadhyay, "A Topological Distance Between Multi-Fields Based on Multi-Dimensional Persistence Diagrams," in IEEE Transactions on Visualization and Computer Graphics, vol. 30, no. 9, pp. 5939-5952, Sept. 2024
