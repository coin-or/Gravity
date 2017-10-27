cmake_minimum_required(VERSION 3.2)

project(Gravity)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Find Gurobi (optional)
option(ENABLE_GUROBI "Link to Gurobi libraries" OFF)

# Find IPOPT (optional)
option(ENABLE_IPOPT "Link to IPOPT libraries" OFF)

# Find BONMIN (optional)
option(ENABLE_BONMIN "Link to BONMIN libraries" OFF)

option(ENABLE_CPLEX "Link to Cplex libraries" ON)

# Find SDPA (optional)
option(ENABLE_SDPA "Link to SDPA libraries" OFF)

# Find MOSEK (optional)
option(ENABLE_MOSEK "Link to MOSEK libraries" OFF)

# include BOOST (optional) 
option(ENABLE_BOOST "Link to BOOST libraries" ON)

# configure a header file to pass some of the CMake settings
# to the source code.
configure_file (
    "${PROJECT_SOURCE_DIR}/GravityConfig.h.in"
    "${PROJECT_BINARY_DIR}/include/gravity/GravityConfig.h")


# Find Ipopt.
if (ENABLE_IPOPT)
	message(STATUS "Enable IPOPT")
	add_definitions(-DUSE_IPOPT)
	find_package(IPOPT REQUIRED)
	include_directories(${IPOPT_INCLUDE_DIRS}/..)
	set(LIBS ${LIBS} ${IPOPT_LIBRARIES})
endif()

# Find SDPA.
if (ENABLE_SDPA)
	message(STATUS "Enable SDPA")
	add_definitions(-DUSE_SDPA)
	find_package(SDPA REQUIRED)
	include_directories(${SDPA_INCLUDE_DIRS})
	set(LIBS ${LIBS} ${SDPA_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-all-loops -Wall -fPIC -funroll-all-loops")
    set(CMAKE_Fortran_FLAGS "-g -O2 -Wall -fPIC -funroll-all-loops")
endif()

# Find MOSEK.
if (ENABLE_MOSEK)
	message(STATUS "Enable MOSEK")
	add_definitions(-DUSE_MOSEK)
	find_package(MOSEK REQUIRED)
	include_directories(${MOSEK_INCLUDE_DIRS})
	set(LIBS ${LIBS} ${MOSEK_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

# Find GUROBI?
if (ENABLE_GUROBI)
    message(STATUS "Enable GUROBI")
    add_definitions(-DUSE_GUROBI)
    find_package(GUROBI)
    include_directories(${GUROBI_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${GUROBI_LIBRARIES})
endif()

# Find BONMIN?
if (ENABLE_BONMIN)
    message(STATUS "Enable BONMIN")
    add_definitions(-DUSE_BONMIN)
    find_package(BONMIN)
    include_directories(${BONMIN_INCLUDE_DIRS}/..)
    set(LIBS ${LIBS} ${BONMIN_LIBRARIES})
endif()

# Find CPLEX?
if (ENABLE_CPLEX)
    message(STATUS "Enable CPLEX")
    add_definitions(-DUSE_CPLEX)
    find_package(CPLEX)
    include_directories(${CPLEX_INCLUDE_DIRS})
    #link_directories(${CPLEX_LINKS})
    set(LIBS ${LIBS} ${CPLEX_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DNDEBUG -DIL_STD")
endif()


# Find BOOST?
if (ENABLE_BOOST)
    message(STATUS "Enable BOOST")
    add_definitions(-DUSE_BOOST)
    include_directories("/usr/local/include/")
    set(LIBS ${LIBS} "usr/local/lib/")
endif()

# Compiler options
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")

# Build options
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

set(CMAKE_BUILD_TYPE "Release")
#set(CMAKE_BUILD_TYPE "Debug")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#set( CMAKE_VERBOSE_MAKEFILE on )

add_subdirectory(src)
add_subdirectory(examples)
# Install
install(FILES ${PROJECT_BINARY_DIR}/PowerToolsConfig.h DESTINATION include/gravity)

# Includes Catch in the project:
#add_subdirectory(src/external/catch)
#include_directories(${CATCH_INCLUDE_DIR} ${COMMON_INCLUDES})
#enable_testing(true)  # Enables unit-testing.
