CMake_Minimum_Required(VERSION 2.8)

Project(Deformetrica)
Set (DEFORMETRICA_VERSION_MAJOR 2)
Set (DEFORMETRICA_VERSION_MINOR 1)

If (COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
Endif (COMMAND cmake_policy)

# == Check gcc compiler: == 
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

enable_testing()

#
# Find packages (Armadillo, CUDA, ITK & VTK) :
#
# Armadillo :
# https://cmake.org/cmake/help/v3.0/module/FindArmadillo.html
Find_Package(Armadillo REQUIRED)
Include_Directories(${ARMADILLO_INCLUDE_DIRS})
IF (ARMADILLO_FOUND)
  Include_directories(${ARMADILLO_INCLUDE_DIRS})
  Message("Armadillo is found")
ELSE (ARMADILLO_FOUND)
  MESSAGE(FATAL_ERROR "Armadillo was not found :(")
ENDIF (ARMADILLO_FOUND)


# ITK :
Find_Package(ITK)
If (ITK_FOUND)
  Include (${ITK_USE_FILE})
Else (ITK_FOUND)
  Message (FATAL_ERROR "ITK was not found :( Please set ITK_DIR.")
Endif (ITK_FOUND)

# VTK :
Find_Package(VTK)
If (VTK_FOUND)
  If ( "${VTK_MAJOR_VERSION}" LESS 6 )
    Message ("VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} is used.")
    Include (${USE_VTK_FILE})
  Else ( "${VTK_MAJOR_VERSION}" LESS 6 )
    Find_Package(VTK COMPONENTS vtkCommonCore vtkCommonDataModel vtkWrappingTools vtkIOLegacy vtkFiltersCore vtkFiltersGeneral NO_MODULE)
    If (VTK_FOUND)
      Message ("VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} is used.")
      Include (${VTK_USE_FILE})
    Else (VTK_FOUND)
      Message (FATAL_ERROR "Some components of VTK have not been found.")
      Endif (VTK_FOUND)
   Endif ( "${VTK_MAJOR_VERSION}" LESS 6 )
Else (VTK_FOUND)
  Message (FATAL_ERROR "VTK was not found :( Please set VTK_DIR.")
Endif (VTK_FOUND)

# Cuda :
Find_Package(CUDA)
If (CUDA_FOUND)
  Option (USE_CUDA "Use cuda implementation of convolutions" ON) 
Else (CUDA_FOUND)
  Message ("CUDA was not found.")
Endif (CUDA_FOUND)

# Release version by default :
If (NOT DEFINED CMAKE_BUILD_TYPE:STRING)
  Set (CMAKE_BUILD_TYPE:STRING Release)
Endif (NOT DEFINED CMAKE_BUILD_TYPE:STRING)

# Basic flags for the C++ compiler :
#Set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -llapack")
#Set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -larmadillo")

# Extra flags for the C++ compiler :
If (CMAKE_BUILD_TYPE MATCHES Release)
  Set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DARMA_NO_DEBUG -DNDEBUG -msse4.2")
Endif (CMAKE_BUILD_TYPE MATCHES Release)

# Extra flags for the C++ compiler if version is Gprof :
If (CMAKE_BUILD_TYPE MATCHES Gprof)
  Set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg")
Endif (CMAKE_BUILD_TYPE MATCHES Gprof)

# Set computations on double precision by default :
If (NOT DEFINED SET_DOUBLE_PRECISION)
  Set (SET_DOUBLE_PRECISION YES)
Endif (NOT DEFINED SET_DOUBLE_PRECISION)
If (SET_DOUBLE_PRECISION)
  Option (USE_DOUBLE_PRECISION "Deformetrica will make computations in double precision" ON) 
Else (SET_DOUBLE_PRECISION)
  Message ("Computations will be in simple precision")
Endif (SET_DOUBLE_PRECISION)

# == Setup Include Directories ==
Include_Directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/deformetrica/src/algo
  ${CMAKE_CURRENT_SOURCE_DIR}/deformetrica/src/algo/atlases
  ${CMAKE_CURRENT_SOURCE_DIR}/deformetrica/src/algo/4D_atlases
  ${CMAKE_CURRENT_SOURCE_DIR}/deformetrica/src/algo/deformable_objects
  ${CMAKE_CURRENT_SOURCE_DIR}/deformetrica/src/algo/deformations
  ${CMAKE_CURRENT_SOURCE_DIR}/deformetrica/src/algo/energy_minimization
  ${CMAKE_CURRENT_SOURCE_DIR}/deformetrica/src/algo/kernels
  ${CMAKE_CURRENT_SOURCE_DIR}/deformetrica/src/algo/linear_algebra
  ${CMAKE_CURRENT_SOURCE_DIR}/deformetrica/src/algo/utilities
  ${CMAKE_CURRENT_SOURCE_DIR}/deformetrica/src/io
  ${CMAKE_CURRENT_SOURCE_DIR}/shape-lme/include
  ${PROJECT_BINARY_DIR}
)

# == What to Build ==
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/deformetrica/)
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/shape-lme/)
