# Configure a header file to pass some of the CMake settings to the source code :
Configure_File(
  "${CMAKE_CURRENT_SOURCE_DIR}/DeformetricaConfig.h.in"
  "${PROJECT_BINARY_DIR}/DeformetricaConfig.h"
)
If (USE_CUDA)

    Message("Deformetrica will be compiled with Cuda")

#    # Set default compute capability to 1.3 in order to enable double precision with Cuda :
#    If (NOT DEFINED COMPUTE_CAPABILITY_GPU)
#        Set (COMPUTE_CAPABILITY_GPU 13)
#    Endif (NOT DEFINED COMPUTE_CAPABILITY_GPU)

#    Set (CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode arch=compute_${COMPUTE_CAPABILITY_GPU},code=sm_${COMPUTE_CAPABILITY_GPU}")

#    Set( MY_CUDA_FLAGS -gencode arch=compute_11,code=sm_11;-gencode arch=compute_12,code=sm_12;-gencode arch=compute_13,code=sm_13;-gencode arch=compute_20,code=sm_20;-gencode arch=compute_30,code=sm_30;-gencode arch=compute_35,code=sm_35)
    Set( MY_CUDA_FLAGS -gencode arch=compute_20,code=sm_20;-gencode arch=compute_30,code=sm_30;-gencode arch=compute_35,code=sm_35)
    List(APPEND CUDA_NVCC_FLAGS ${MY_CUDA_FLAGS})

    # Library containing the Cuda operations
    Set (LIB_TYPE STATIC)
    Cuda_Add_Library (
    CudaConvolution ${LIB_TYPE}
        ${CMAKE_CURRENT_SOURCE_DIR}/../lib/cuda_convolutions/RadialFunction.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../lib/cuda_convolutions/ScalarRadialKernel.h
        ${CMAKE_CURRENT_SOURCE_DIR}/../lib/cuda_convolutions/GpuConv1D.cu
    )
    Set (CUDA_LIBS ${CUDA_LIBS} CudaConvolution)

Cuda_Add_Executable(benchmark
  ../src/algo/utilities/SimpleTimer.cxx
  cuda_benchmark/benchmark.cxx
)

Cuda_Add_Executable(sparseMatching2
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  matching/sparseMatching2.cxx
)

Cuda_Add_Executable(sparseMatching3
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  matching/sparseMatching3.cxx
)

Cuda_Add_Executable(sparseAtlas2
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  atlas_building/sparseAtlas2.cxx
)

Cuda_Add_Executable(sparseAtlas3
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  atlas_building/sparseAtlas3.cxx
)

Cuda_Add_Executable(sparseGeodesicRegression2
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  regression/sparseRegression2.cxx
)

Cuda_Add_Executable(sparseGeodesicRegression3
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  regression/sparseRegression3.cxx
)

Cuda_Add_Executable(ShootAndFlow2
	../src/io/SparseDiffeoParametersXMLFile.cxx
	../src/io/DeformableObjectParametersXMLFile.cxx
	../src/io/SparseDiffeoParameters.cxx
	../src/io/DeformableObjectParameters.cxx
	../src/algo/utilities/myvtkPolyDataNormals.cxx
	shooting/ShootAndFlow2.cxx
)

Cuda_Add_Executable(ShootAndFlow3
	../src/io/SparseDiffeoParametersXMLFile.cxx
	../src/io/DeformableObjectParametersXMLFile.cxx
	../src/io/SparseDiffeoParameters.cxx
	../src/io/DeformableObjectParameters.cxx
	../src/algo/utilities/myvtkPolyDataNormals.cxx
	shooting/ShootAndFlow3.cxx
)

Cuda_Add_Executable(CreateDeformationField
    ../src/io/DeformationFieldIO.txx
	../src/io/SparseDiffeoParameters.cxx
	../src/io/SparseDiffeoParametersXMLFile.cxx
	../src/io/DeformableObjectParameters.cxx
	../src/algo/utilities/myvtkPolyDataNormals.cxx
	deformation/CreateDeformationField.cxx
)


Target_Link_Libraries(benchmark
  ${ARMADILLO_LIBRARIES}
  ${ITK_LIBRARIES}
  ${CUDA_LIBS}
)



else(USE_CUDA)

    Message("Deformetrica will be compiled without Cuda")

Add_Executable(sparseMatching2
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  matching/sparseMatching2.cxx
)

Add_Executable(sparseMatching3
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  matching/sparseMatching3.cxx
)

Add_Executable(sparseAtlas2
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  atlas_building/sparseAtlas2.cxx
)

Add_Executable(sparseAtlas3
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  atlas_building/sparseAtlas3.cxx
)

Add_Executable(sparseGeodesicRegression2
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  regression/sparseRegression2.cxx
)

Add_Executable(sparseGeodesicRegression3
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/SimpleTimer.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  regression/sparseRegression3.cxx
)

Add_Executable(ShootAndFlow2
    ../src/io/SparseDiffeoParametersXMLFile.cxx
    ../src/io/DeformableObjectParametersXMLFile.cxx
    ../src/io/SparseDiffeoParameters.cxx
    ../src/io/DeformableObjectParameters.cxx
    ../src/algo/utilities/myvtkPolyDataNormals.cxx
    shooting/ShootAndFlow2.cxx
)

Add_Executable(ShootAndFlow3
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParametersXMLFile.cxx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  shooting/ShootAndFlow3.cxx
)

Add_Executable(CreateDeformationField
  ../src/io/DeformationFieldIO.txx
  ../src/io/SparseDiffeoParameters.cxx
  ../src/io/SparseDiffeoParametersXMLFile.cxx
  ../src/io/DeformableObjectParameters.cxx
  ../src/algo/utilities/myvtkPolyDataNormals.cxx
  deformation/CreateDeformationField.cxx
)

endif(USE_CUDA)


# Add_Executable(CreateDeformationField
#   ../src/io/SparseDiffeoParametersXMLFile.cxx
#   ../src/io/DeformationFieldIO.txx
#   ../src/io/SparseDiffeoParameters.cxx
#   ../src/algo/utilities/myvtkPolyDataNormals.cxx
#   matching/CreateDeformationField.cxx
# )



#
# Target link :
#
Set ( VTK_DEPENDANCIES "" )
If ( "${VTK_MAJOR_VERSION}" LESS 6 )
    List (APPEND VTK_DEPENDANCIES
        vtkCommon
        vtkGraphics
        vtkIO
        vtkRendering
        vtkFiltering
    )
Else ( "${VTK_MAJOR_VERSION}" LESS 6 )
    List (APPEND VTK_DEPENDANCIES
        ${VTK_LIBRARIES}
    )
Endif ( "${VTK_MAJOR_VERSION}" LESS 6 )


Target_Link_Libraries(sparseMatching2
  ${ARMADILLO_LIBRARIES}
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(sparseMatching3
  ${ARMADILLO_LIBRARIES}
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(sparseAtlas2
  ${ARMADILLO_LIBRARIES}
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(sparseAtlas3
  ${ARMADILLO_LIBRARIES}
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(sparseGeodesicRegression2
  ${ARMADILLO_LIBRARIES}
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(sparseGeodesicRegression3
  ${ARMADILLO_LIBRARIES}
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)


Target_Link_Libraries(ShootAndFlow2
  ${ARMADILLO_LIBRARIES}
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(ShootAndFlow3
  ${ARMADILLO_LIBRARIES}
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)

Target_Link_Libraries(CreateDeformationField
  ${ARMADILLO_LIBRARIES}
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
  ${CUDA_LIBS}
)



#
# Dependencies (only for VTK < 6) :
#
If ( "${VTK_MAJOR_VERSION}" LESS 6 )
Add_Dependencies (sparseMatching2
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (sparseMatching3
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (sparseAtlas2
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (sparseAtlas3
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (sparseGeodesicRegression2
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (sparseGeodesicRegression3
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (ShootAndFlow2
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (ShootAndFlow3
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Add_Dependencies (CreateDeformationField
   ${ITK_LIBRARIES}
   ${VTK_DEPENDANCIES}
)

Endif ( "${VTK_MAJOR_VERSION}" LESS 6 )





# add tests
enable_testing()

add_test(NAME surface_matching WORKING_DIRECTORY ../../examples/surface_matching/ COMMAND sparseMatching3 paramDiffeos.xml paramSurface.xml sourceSurfaceMesh.vtk targetSurfaceMesh.vtk)





