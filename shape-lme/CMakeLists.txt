PROJECT(ShapeLME)

CMAKE_MINIMUM_REQUIRED(VERSION 2.7)

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

# == Adds any local CMake models (Find*.cmake files): ==
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${ShapeLME_SOURCE_DIR}/CMake/Modules/")

# == Using Armadillo: ==
FIND_PACKAGE(Armadillo REQUIRED)
IF (ARMADILLO_FOUND)
  INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})
  set(ARMADILLO_INCLUDES ${ARMADILLO_INCLUDE_DIRS})
ELSE()
  message("ARMADILLO NOT FOUND")
ENDIF()

# == Optional subdirectories to build ==
OPTION(BUILD_TESTING "Build tests?" ON)

# == Setup Include Directories ==
LIST(APPEND ShapeLME_INCLUDE_DIRECTORIES
  ${ShapeLME_SOURCE_DIR}/include/)
INCLUDE_DIRECTORIES(${ShapeLME_INCLUDE_DIRECTORIES})

SET(ShapeLME
  ${ARMADILLO_LIBRARIES})
LINK_DIRECTORIES(${ShapeLME})

# == Shared Library Option ==
MARK_AS_ADVANCED(BUILD_SHARED_LIBS)
OPTION(BUILD_SHARED_LIBS "Build shared libraries?" ON)
IF(BUILD_SHARED_LIBS)
  SET(LIBTYPE SHARED)
ELSE(BUILD_SHARED_LIBS)
  SET(LIBTYPE STATIC)
ENDIF(BUILD_SHARED_LIBS)

# == What to Build ==
IF(BUILD_TESTING)
ADD_EXECUTABLE(hypothesis-testing test/hypothesis-testing.cxx)
TARGET_LINK_LIBRARIES(hypothesis-testing ${ShapeLME})

ADD_EXECUTABLE(mixed-effects-covariate-estimation test/mixed-effects-covariate-estimation.cxx)
TARGET_LINK_LIBRARIES(mixed-effects-covariate-estimation ${ShapeLME})
ENDIF(BUILD_TESTING)
