set (CMAKE_CXX_STANDARD 11)
project(META)

cmake_minimum_required(VERSION 2.8)

# Find and load VTK settings.
if(NOT VTK_BINARY_DIR)
  find_package(VTK 6.3 REQUIRED NO_MODULE)
  include(${VTK_USE_FILE})
endif()

if(VTK_WRAP_PYTHON)
  option(BUILD_PYTHON_WRAPPERS "Build python wrappers" ON)
else()
  unset(BUILD_PYTHON_WRAPPERS CACHE)
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Choose static or shared libraries.
option(BUILD_SHARED_LIBS "Build shared libraries." Off)

# ---------------------------------------------------------------------
# IMPORTANT
# You need to provide three paths:
# 1.  The path to the root directory of your OGDF build.
# 2.  The path to your boost library include directory.
# 3.  The path to the subdirectory of your OGDF build containing libOGDF.a
# In each case below, replace the text delimited by <...> (including the
# delimiters) with the path.

include_directories(SYSTEM
  ${VTK_INCLUDE_DIRS}
  ${Boost_LIBRARIES}
  /home/iiitb/Softwares/OGDF/
  /home/iiitb/Softwares/OGDF/include/
)
link_directories(${LINK_DIRECTORIES} /home/iiitb/Softwares/OGDF/_release/ ) 

# ---------------------------------------------------------------------



# Set up our directory structure for output libraries and binaries
# (Note: these are the build locations, not the install locations)
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${META_BINARY_DIR}/bin")
endif()
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  if(UNIX)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${META_BINARY_DIR}/lib")
  else()
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${META_BINARY_DIR}/bin")
  endif()
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${META_BINARY_DIR}/lib")
endif()

# Set up our target directory structure for "make install"
set(META_BIN_DIR "/bin") # for executables and ".dll" libraries
set(META_LIB_DIR "/lib") # for unix/linux/osx shared libraries
set(META_ARC_DIR "/lib") # for static or ".lib" libraries
set(META_INC_DIR "/include") # for header files
set(META_DATA_DIR "/share/meta-${META_SHORT_VERSION}")
set(META_CMAKE_DIR "/lib/meta-${META_SHORT_VERSION}/cmake")


# Set the full paths to the install tree
set(META_RUNTIME_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}${META_BIN_DIR})
set(META_LIBRARY_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}${META_LIB_DIR})
set(META_ARCHIVE_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}${META_ARC_DIR})
set(META_INCLUDE_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}${META_INC_DIR})
set(META_DATA_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}${META_DATA_DIR})
set(META_CMAKE_INSTALL_DEST ${CMAKE_INSTALL_PREFIX}${META_CMAKE_DIR})

# Include directories
set(META_INCLUDE_DIRS "${META_SOURCE_DIR}" "${META_BINARY_DIR}")
include_directories(${DICOM_INCLUDE_DIRS})

# Create the main library
set(LIB_NAME vtkMETA)



#if (CMAKE_INSTALL_PREFIX) 
#  message("SET")
#  message(${CMAKE_INSTALL_PREFIX})
#else()
#  message("NOT SET")
#endif()

# Source files.
set(LIB_SRCS
vtkExtractReebGraph.cxx
vtkJointContourNet.cxx
vtkMultiDimensionalReebGraph.cxx
vtkPolytopeGeometry.cxx
vtkPyramidTree.cxx
vtkSimplicate.cxx
vtkReduceGraph.cxx
vtkLocalOverrides.cxx
vtkMultiResolutionReebSpace.cxx
vtkSimilarityBetweenMultiResolutionReebSpaces.cxx
vtkJointContourNetForShapes.cxx
vtkComputeUnstructuredGrid.cxx
vtkPolytopeGeometryForShapes.cxx
vtkHandleDegeneracyInReebGraph.cxx
vtkDistanceBetweenMultiDimensionalReebGraphs.cxx
vtkPersistenceDiagramOfReebGraph.cxx
vtkWassersteinDistanceBetweenMultiDimensionalReebGraphs.cxx
vtkJCNSplitField.cxx
vtkRAWReader.cxx
)
#vtkMDRGSimilarityByPersistenceStandardBottleNeckDistance.cxx
# Automatically find the header for each cxx file
# NOT COPIED

# hierarchy file
if(BUILD_PYTHON_WRAPPERS AND "${VTK_MAJOR_VERSION}" GREATER 5)
  set_source_files_properties(${LIB_SPECIAL} PROPERTIES WRAP_SPECIAL ON)
  set(vtk-module ${LIB_NAME})
  set(MODULE_HIERARCHY_NAME ${LIB_NAME}Hierarchy)
#  set(${LIB_NAME}_LINK_DEPENDS ${VTK_LIBRARIES})
#  include(${VTK_CMAKE_DIR}/vtkWrapHierarchy.cmake)
#  vtk_wrap_hierarchy(${MODULE_HIERARCHY_NAME} ${META_BINARY_DIR} "${LIB_SRCS}")
#  set(KIT_HIERARCHY_FILE ${META_BINARY_DIR}/${MODULE_HIERARCHY_NAME}.txt)
#  set(LIB_HIERARCHY_STAMP ${META_BINARY_DIR}/${MODULE_HIERARCHY_NAME}.stamp)
endif()

# Create the libary
if(VTK_LIBRARIES)
  set(VTK_LIBS ${VTK_LIBRARIES}) # VTK 6.x
else()
  set(VTK_LIBS vtkIO) # VTK 5.x
endif()
add_library(${LIB_NAME} ${LIB_SRCS} ${LIB_HIERARCHY_STAMP})
target_link_libraries(${LIB_NAME} ${VTK_LIBS} OGDF)
set(META_LIBRARIES ${LIB_NAME})

# Wrappers
if(BUILD_PYTHON_WRAPPERS)
  set(vtk-module ${LIB_NAME})
  set(MODULE_PYTHON_NAME ${LIB_NAME}Python)
  set(LIB_PYTHON_NAME ${LIB_NAME}PythonD)
  if("${VTK_MAJOR_VERSION}" GREATER 5)
    set(LIB_PYTHON_LIBS vtkCommonCorePythonD vtkIOImagePythonD vtkIOLegacyPythonD vtkInfovisCorePythonD vtkInfovisLayoutPythonD vtkInteractionStylePythonD)
  else()
    set(LIB_PYTHON_LIBS vtkCommonPythonD vtkIOPythonD)
  endif()

  if(VTK_PYTHON_INCLUDE_DIR)
    include_directories("${VTK_PYTHON_INCLUDE_DIR}")
  endif()
  include(${VTK_CMAKE_DIR}/vtkWrapPython.cmake)
  vtk_wrap_python3(${MODULE_PYTHON_NAME} LIB_PYTHON_SRCS "${LIB_SRCS}")
  add_library(${LIB_PYTHON_NAME} ${LIB_PYTHON_SRCS} ${LIB_PYTHON_EXTRA_SRCS})
  target_link_libraries(${LIB_PYTHON_NAME} ${LIB_NAME} ${LIB_PYTHON_LIBS} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
  #
  # On Win32 and Mac, link python library now, on UNIX, don't
  if(WIN32 OR APPLE)
    target_link_libraries(${LIB_PYTHON_NAME} ${VTK_PYTHON_LIBRARIES} ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
  endif()
  if(KIT_PYTHON_DEPS)
    add_dependencies(${LIB_PYTHON_NAME} ${KIT_PYTHON_DEPS})
  endif()
  add_library(${MODULE_PYTHON_NAME} MODULE ${MODULE_PYTHON_NAME}Init.cxx)
  set_target_properties(${MODULE_PYTHON_NAME} PROPERTIES PREFIX "")
  if(WIN32 AND NOT CYGWIN)
    set_target_properties(${MODULE_PYTHON_NAME} PROPERTIES SUFFIX ".pyd")
  endif()
  target_link_libraries(${MODULE_PYTHON_NAME} ${LIB_PYTHON_NAME})
  # On Win32 and Mac, link python library now, on UNIX, don't
  if(WIN32 OR APPLE)
    target_link_libraries(${MODULE_PYTHON_NAME} ${VTK_PYTHON_LIBRARIES})
  endif()
endif()

include_directories(${META_INCLUDE_DIRS} ${META_INCLUDE_DIRS}/override)

# Set the install rules for the library
install(TARGETS ${LIB_NAME} ${LIB_PYTHON_NAME} ${MODULE_PYTHON_NAME}
  EXPORT METATargets
  RUNTIME DESTINATION ${META_RUNTIME_INSTALL_DEST} COMPONENT RuntimeLibraries
  LIBRARY DESTINATION ${META_LIBRARY_INSTALL_DEST} COMPONENT RuntimeLibraries
  ARCHIVE DESTINATION ${META_ARCHIVE_INSTALL_DEST} COMPONENT Development)

install(FILES ${LIB_HDRS}
  DESTINATION ${META_INCLUDE_INSTALL_DEST} COMPONENT Development)



