cmake_minimum_required(VERSION 3.2)

project(Gravity)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
set(MACOSX_RPATH 1) 
# Debug mode
option(Debug "Turn on debug mode" OFF)

# Coverage
option(Cov "Turn on coverage" OFF)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(Matplot "Link to Matplotlib libraries" OFF)
if (Matplot)
 add_definitions(-DUSE_MATPLOT)
 find_package(Python3 COMPONENTS Development NumPy)
endif()

# Find GUROBI (optional)
option(Gurobi "Link to GUROBI libraries" OFF)

# Find IPOPT (REQUIRED)
option(Ipopt "Link to IPOPT libraries" ON)

# Find LASlib (REQUIRED)
option(LASlib "Link to LASlib libraries" ON)

# Find QHULL (OPTIONAL)
option(QHULL "Link to QHULL libraries" OFF)

# Find VORO (REQUIRED)
option(VORO "Link to VORO libraries" ON)

# Find XLNT (optional)
option(Xlnt "Link to XLNT libraries" OFF)

# Find BONMIN (optional)
option(Bonmin "Link to BONMIN libraries" OFF)

# Find CPLEX (optional)
option(Cplex "Link to CPLEX libraries" OFF)

# Find SDPA (optional)
option(Sdpa "Link to SDPA libraries" OFF)

# Find MOSEK (optional)
option(Mosek "Link to MOSEK libraries" OFF)

# include BOOST (optional)
option(Boost "Link to BOOST libraries" OFF)

# include QPP (optional)
option(Qpp "Link to QPP libraries" OFF)

# Find CLP (optional)
option(Clp "Link to CLP libraries" OFF)

# Find OpenMPI (optional)
option(OpenMPI "Link to OpenMPI libraries" OFF)

# Find OPT_PARSER (optional)
option(OPT_PARSER "Link to OPT_PARSER libraries" OFF)

if (OpenMPI)
    SET(CMAKE_CXX_COMPILER mpicxx)
    SET(CMAKE_C_COMPILER mpicc)
    message("Compilers used: ${CMAKE_C_COMPILER} and ${CMAKE_CXX_COMPILER}")
    message(STATUS "Enable OpenMPI")
    add_definitions(-DUSE_MPI)
endif()



# Find Ipopt.
if (Ipopt)
message(STATUS "Enable IPOPT")
add_definitions(-DUSE_IPOPT)
find_package(IPOPT QUIET)
if (IPOPT_FOUND)
 message("-- Found IPOPT: ${IPOPT_INCLUDE_DIRS}")
 include_directories(${IPOPT_INCLUDE_DIRS})
 set(LIBS ${LIBS} ${IPOPT_LIBRARIES})
else(IPOPT_FOUND)
set(ADD_IPOPT TRUE)
if(APPLE)
 set(IPOPT_LIBRARIES ${THIRDPARTY_INSTALL_PATH}/Install/ipopt/build/lib/libipopt.dylib)
elseif(UNIX)
 set(IPOPT_LIBRARIES ${THIRDPARTY_INSTALL_PATH}/Install/ipopt/build/lib/libipopt.so)
endif(APPLE)
 include(ExternalIpopt)
endif (IPOPT_FOUND)
endif()
# Find xlnt
if (Xlnt)
	message(STATUS "Enable XLNT")
	add_definitions(-DUSE_XLNT)
	find_package(XLNT REQUIRED)
	include_directories(${XLNT_INCLUDE_DIRS}/..)
	set(LIBS ${LIBS} ${XLNT_LIBRARIES})
endif()


# Find SDPA.
if (Sdpa)
	message(STATUS "Enable SDPA")
	add_definitions(-DUSE_SDPA)
	find_package(SDPA REQUIRED)
	include_directories(${SDPA_INCLUDE_DIRS})
	set(LIBS ${LIBS} ${SDPA_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-all-loops -Wall -fPIC -funroll-all-loops")
    set(CMAKE_Fortran_FLAGS "-g -O2 -Wall -fPIC -funroll-all-loops")
endif()

# Find MOSEK.
if (Mosek)
    set(MOSEK_ROOT_DIR "$ENV{MOSEK_ROOT_DIR}" CACHE PATH "MOSEK root directory.")
	message(STATUS "Enable MOSEK")
	add_definitions(-DUSE_MOSEK)
	find_package(MOSEK REQUIRED)
	include_directories(${MOSEK_INCLUDE_DIRS})
	set(LIBS ${LIBS} ${MOSEK_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    link_directories(${MOSEK_ROOT_DIR}/bin)
    link_directories(${MOSEK_ROOT_DIR}/src/fusion_cxx)

endif()

# Find gurobi?
if (Gurobi)
    message(STATUS "Enable GUROBI")
    add_definitions(-DUSE_GUROBI)
    find_package(GUROBI)
    include_directories(${GUROBI_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${GUROBI_LIBRARIES})
endif()

# Find CLP?
if (Clp)
    message(STATUS "Enable CLP")
    add_definitions(-DUSE_CLP)
    find_package(CLP REQUIRED)
    include_directories(${CLP_INCLUDE_DIRS}/)
    set(LIBS ${LIBS} ${CLP_LIBRARIES})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lCoinUtils -lClp")
endif()

# Find BONMIN?
if (Bonmin)
    message(STATUS "Enable BONMIN")
    add_definitions(-DUSE_BONMIN)
    find_package(BONMIN)
    include_directories(${BONMIN_INCLUDE_DIRS})
    include_directories(${BONMIN_INCLUDE_DIRS}/..)
    set(LIBS ${LIBS} ${BONMIN_LIBRARIES})
endif()

# Find CPLEX?
if (Cplex)
    message(STATUS "Enable CPLEX")
    add_definitions(-DUSE_CPLEX)
    find_package(CPLEX)
    include_directories(${CPLEX_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${CPLEX_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DNDEBUG -DIL_STD")
endif()


# Find BOOST?
if (Boost)
    message(STATUS "Enable BOOST")
    add_definitions(-DUSE_BOOST)
    include_directories("/usr/local/include/")
    set(LIBS ${LIBS} "usr/local/lib/")
endif()


set(UNZIPPED OFF CACHE BOOL "Already unzipped")
if(NOT UNZIPPED)
  message(STATUS "Unzipping Classification Data")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/data_sets/Classification/Archive.tar.gz
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data_sets/Classification
     )
endif()
set_directory_properties(PROPERTIES EP_BASE ${CMAKE_BINARY_DIR}/third_party)
get_directory_property(THIRDPARTY_BASE_PATH EP_BASE)

set(THIRDPARTY_INSTALL_PATH ${THIRDPARTY_BASE_PATH})


include(ExternalProject)
# Eigen
#find_package(Eigen QUIET)
#if (EIGEN_FOUND)
# message("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
#else(EIGEN_FOUND)
# set(EIGEN_VERSION 3.3.7)
# include(ExternalEigen)
# set(EIGEN_INCLUDE_DIR ${THIRDPARTY_INSTALL_PATH}/Install/eigen)
# set(ADD_EIGEN TRUE)
#endif (EIGEN_FOUND)
#include_directories(${EIGEN_INCLUDE_DIR})

#find_package(GoICP QUIET)
#if (GoICP_FOUND)
# message("-- Found GoICP: ${EIGEN_INCLUDE_DIRS}")
#else(GoICP_FOUND)
# include(ExternalGoICP)
# set(GoICP_INCLUDE_DIR ${THIRDPARTY_INSTALL_PATH}/Install/GoICP)
# set(ADD_GoICP TRUE)
#endif (GoICP_FOUND)
#include_directories(${GoICP_INCLUDE_DIR})

if (QHULL)
message(STATUS "Enable QHULL")
add_definitions(-DUSE_QHULL)
find_package(QHULL QUIET)
if (QHULL_FOUND)
 message("-- Found QHULL: ${QHULL_INCLUDE_DIRS}")
 include_directories(${QHULL_INCLUDE_DIRS})
 set(LIBS ${LIBS} ${QHULL_LIBRARIES})
else(QHULL_FOUND)
 set(ADD_QHULL TRUE)
 message(STATUS "Will download QHULL library from http://www.qhull.org/download")
 include(ExternalQHULL)
 #set(QHULL_LIBRARIES ${THIRDPARTY_INSTALL_PATH}/Install/QHULL/build/libqhullcpp.a ${THIRDPARTY_INSTALL_PATH}/Install/QHULL/build/libqhullstatic.a ${THIRDPARTY_INSTALL_PATH}/Install/QHULL/build/libqhullstatic_r.a)
 set(QHULL_LIBRARIES ${THIRDPARTY_INSTALL_PATH}/Install/qhull/build/libqhullcpp.a ${THIRDPARTY_INSTALL_PATH}/Install/qhull/build/libqhullstatic_r.a)
endif()
endif()

if (VORO)
message(STATUS "Enable VORO")
add_definitions(-DUSE_VORO)
find_package(VORO QUIET)
if (VORO_FOUND)
 message("-- Found VORO: ${VORO_INCLUDE_DIRS}")
 include_directories(${VORO_INCLUDE_DIRS})
 set(LIBS ${LIBS} ${VORO_LIBRARIES})
else(VORO_FOUND)
 set(ADD_VORO TRUE)
 message(STATUS "Will download VORO library from http://math.lbl.gov/voro++/download/")
 include(ExternalVORO)
 set(VORO_LIBRARIES ${THIRDPARTY_INSTALL_PATH}/Install/voro/src/libvoro++.a) 
endif()
endif()

if (LASlib)
message(STATUS "Enable LASlib")
add_definitions(-DUSE_LASlib)
find_package(LASlib QUIET)
if (LASlib_FOUND)
 message("-- Found LASlib: ${LASlib_INCLUDE_DIRS}")
 include_directories(${LASlib_INCLUDE_DIRS})
 set(LIBS ${LIBS} ${LASlib_LIBRARIES})
else(LASlib_FOUND)
 set(ADD_LASlib TRUE)
 message(STATUS "Will download LASlib library from https://github.com/LAStools/LAStools")
 include(ExternalLASTOOLS)
 set(LASlib_LIBRARIES ${THIRDPARTY_INSTALL_PATH}/Install/lastools/LASlib/lib/liblas.a)
endif()
endif()

# cpp_option_parser
if (OPT_PARSER)
	set(OPT_PARSER_VERSION v1.0)
	message(STATUS "Enable OPT_PARSER")
	add_definitions(-DUSE_OPT_PARSER)
	include(ExternalOptParser)
	set(OPT_PARSER_INCLUDE_DIR ${THIRDPARTY_INSTALL_PATH}/Source/opt_parser)
	link_directories(${THIRDPARTY_INSTALL_PATH}/Build/opt_parser)
	link_directories(${THIRDPARTY_INSTALL_PATH}/Build/opt_parser/Debug)
	include_directories(${OPT_PARSER_INCLUDE_DIR})
	include_directories(${THIRDPARTY_INSTALL_PATH}/include)
endif()



# Find QPP?
if (Qpp AND NOT UNZIPPED)
message(STATUS "Unzipping QPP v1.0")
execute_process(
  COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/third_party/qpp-1.0.tar.gz
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/third_party
  )
    add_definitions(-DUSE_QPP)
    include_directories(${PROJECT_SOURCE_DIR}/third_party/qpp-1.0/include)
endif()

set(UNZIPPED ON CACHE BOOL "Already unzipped" FORCE)
# Inlude dir
include_directories(${PROJECT_SOURCE_DIR}/third_party)


# Include ARMO dirs
include_directories(${PROJECT_SOURCE_DIR}/ARMO)

set(LIBS ${LIBS} ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/lib/Release ${PROJECT_SOURCE_DIR}/lib/Debug)
# Compiler options
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")

# Build options
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(PROJECT_BINARY_DIR ${PROJECT_SOURCE_DIR}/build)

set(CMAKE_BUILD_TYPE "Release")
if (Debug)
set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
endif()

if(Cov)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -fprofile-arcs -lgcov")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'\$ORIGIN/lib'")


if(LINUX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()


add_subdirectory(src)
add_subdirectory(ARMO)

#Get user defined configuration including passing project source directory into c++ code
configure_file(GravityConfig.h.in ${PROJECT_SOURCE_DIR}/include/gravity/GravityConfig.h)
