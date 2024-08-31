# Start with an Ubuntu 16.04 base image
FROM ubuntu:16.04

# Install necessary dependencies and update cmake
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y \
    python2.7 \
    build-essential \
    libboost-all-dev \
    libpython2.7-dev \
    software-properties-common \
    wget \
    libgl1-mesa-dev \
    vim \
    libglu1-mesa-dev \
    freeglut3-dev \
    libglew-dev && \
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | apt-key add - && \
    apt-add-repository 'deb https://apt.kitware.com/ubuntu/ xenial main' && \
    apt-get update && \
    apt-get install -y cmake && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy VTK source file and OGDF folder
COPY VTK-6.3.0 /app/VTK-6.3.0
COPY OGDF /app/ogdf




RUN cd VTK-6.3.0 && mkdir build && cd build && cmake -DBUILD_SHARED_LIBS=ON     -DVTK_WRAP_PYTHON=ON   -DModule_vtkFiltersReebGraph=ON -DModule_vtkInfovisBoost=ON -DModule_vtkInfovisBoostGraphAlgorithms=ON  -DModule_vtkPythonInterpreter=ON  .. && make -j10 && make install

# Build OGDF
RUN chmod +x /app/ogdf/makeMakefile.sh
RUN cd /app/ogdf && \
    sed -i 's/sharedLib false/sharedLib true/' makeMakefile.config&& \
    sed -i 's/libCoin true/libCoin false/' makeMakefile.config && \
    ./makeMakefile.sh && \
    make -j20 && make install
COPY vtkLocal  /app/vtkLocal
COPY Python /app/Python

# Build vtkLocal
RUN cd  vtkLocal && \
    mkdir build && cd build && \
    cmake -DCMAKE_CXX_FLAGS="-I/app/ogdf/include/ogdf" -DCMAKE_EXE_LINKER_FLAGS="-L/app/ogdf/_release -lOGDF"  -DCMAKE_MODULE_LINKER_FLAGS="-L/app/ogdf/_release -lOGDF" -DCMAKE_SHARED_LINKER_FLAGS="-L/app/ogdf/_release -lOGDF"  .. && \
    make -j20 && make install

# Set environment variables
ENV LOCALPATH=/app/vtkLocal
ENV VTKPATH=/app/VTK-6.3.0/build
ENV VTKCMAKEPATH=/usr/local/lib/cmake/vtk-6.3
ENV PYTHONPATH=$LOCALPATH:$LOCALPATH/lib:/usr/local/include/vtk-6.3:$VTKPATH:$VTKPATH/Wrapping/Python:/usr/local/lib/python2.7/site-packages/vtk:/usr/local/lib
ENV LD_LIBRARY_PATH=/app/VTK-6.3.0/lib:/app/vtkLocal/lib:/usr/local/lib:$LOCALPATH:$LOCALPATH/lib:$VTKPATH/bin:$VTKCMAKEPATH:/usr/local/lib/python2.7/site-packages/vtk
ENV PATH=/usr/local/lib/python2.7/site-packages/vtk:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/local/lib


#Copy the test data and the shell script file for computing the distances
COPY TestData /app/data/
COPY compute-distance.sh /app/compute-distance.sh
RUN chmod +x /app/compute-distance.sh

ENTRYPOINT ["/app/compute-distance.sh"]


# Run this for more customizations
# ENTRYPOINT ["/bin/bash"]   
