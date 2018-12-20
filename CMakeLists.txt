cmake_minimum_required(VERSION 3.2)

project(Gravity)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Debug mode
option(debug "Turn on debug mode" OFF)

# Find GUROBI (optional)
option(Gurobi "Link to GUROBI libraries" OFF)

# Find IPOPT (optional)
option(Ipopt "Link to IPOPT libraries" OFF)

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


# Find Ipopt.
if (Ipopt)
	message(STATUS "Enable IPOPT")
	add_definitions(-DUSE_IPOPT)
	find_package(IPOPT REQUIRED)
	include_directories(${IPOPT_INCLUDE_DIRS}/..)
	set(LIBS ${LIBS} ${IPOPT_LIBRARIES})
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
	message(STATUS "Enable MOSEK")
	add_definitions(-DUSE_MOSEK)
	find_package(MOSEK REQUIRED)
	include_directories(${MOSEK_INCLUDE_DIRS})
	set(LIBS ${LIBS} ${MOSEK_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
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
# Include Eigen

set(UNZIPPED OFF CACHE BOOL "Already unzipped")
if(NOT UNZIPPED)
  message(STATUS "Unzipping Classification Data and Eigen v3.3.5")
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/data_sets/Classification/Archive.tar.gz
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data_sets/Classification
     )
execute_process(
  COMMAND ${CMAKE_COMMAND} -E tar xzf ${PROJECT_SOURCE_DIR}/third_party/eigen-git-mirror-3.3.5.tar.gz
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/third_party
  )
endif()
include_directories(${PROJECT_SOURCE_DIR}/third_party/eigen-git-mirror-3.3.5)

# Include cpp_option_parser
include_directories(${PROJECT_SOURCE_DIR}/third_party/cpp_option_parser)
add_subdirectory(third_party/cpp_option_parser)
set(LIBS ${LIBS} option_parser)

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

# Include Power dirs
include_directories(${PROJECT_SOURCE_DIR}/examples/Optimization/NonLinear/Power)

# Include ML dirs
include_directories(${PROJECT_SOURCE_DIR}/examples/MachineLearning)

set(LIBS ${LIBS} ${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/lib/Release ${PROJECT_SOURCE_DIR}/lib/Debug)
# Compiler options
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")

# Build options
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Release)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

set(CMAKE_BUILD_TYPE "Release")
if (debug)
set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin/Debug)
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lpthread  -fprofile-arcs -ftest-coverage")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'\$ORIGIN/lib'" -fprofile-arcs -lgcov)

add_subdirectory(src)
add_subdirectory(examples)

#Get user defined configuration including passing project source directory into c++ code
configure_file(GravityConfig.h.in ${PROJECT_SOURCE_DIR}/include/gravity/GravityConfig.h)
