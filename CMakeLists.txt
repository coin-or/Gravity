cmake_minimum_required(VERSION 3.2)

project(Gravity)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Find solvers on system
#find_package(GUROBI REQUIRED)
#find_package(IPOPT REQUIRED)
#find_package(BONMIN REQUIRED)

# Include solver headers and libraries
#include_directories(${GUROBI_INCLUDE_DIRS} ${IPOPT_INCLUDE_DIRS}/.. ${BONMIN_INCLUDE_DIRS}/..)
#set(LIBS ${LIBS} ${GUROBI_LIBRARIES} ${IPOPT_LIBRARIES} ${BONMIN_LIBRARIES})

# Compiler options

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast")

# Build options
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

#set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_BUILD_TYPE "Debug")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#set( CMAKE_VERBOSE_MAKEFILE on )

add_subdirectory(src)

# Includes Catch in the project:
add_subdirectory(src/external/catch)
include_directories(${CATCH_INCLUDE_DIR} ${COMMON_INCLUDES})
enable_testing(true)  # Enables unit-testing.
