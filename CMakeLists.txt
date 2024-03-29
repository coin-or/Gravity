cmake_minimum_required(VERSION 3.16)

project(Gravity)

enable_testing()

include(ExternalProject)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
set(MACOSX_RPATH 1) 
# Debug mode
option(Debug "Turn on debug mode" OFF)

# Coverage
option(Cov "Turn on coverage" OFF)

# Find PcapPlusPlus (optional)
option(PcapPlusPlus "Link to PcapPlusPlus libraries" OFF)

# Find GUROBI (optional)
option(Gurobi "Link to GUROBI libraries" OFF)

# Find HiGHS (optional)
option(HiGHS "Link to HiGHS libraries" OFF)

# Find IPOPT (optional)
option(Ipopt "Link to IPOPT libraries" ON)

option(MP "Link to AMPL MP libraries" OFF)

option(CoinUtils "Link to CoinUtils libraries" OFF)

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
#   find_package(OpenMPI REQUIRED)
#   include_directories(${OpenMPI_INCLUDE_DIRS})
#   set(LIBS ${LIBS} ${OpenMPI_LIBRARIES})
endif()


if (PcapPlusPlus)
    message(STATUS "Enable PcapPlusPlus")
    add_definitions(-DUSE_PcapPlusPlus)
    find_package(PcapPlusPlus REQUIRED)
    include_directories(${PcapPlusPlus_INCLUDE_DIRS}/pcapplusplus)
    set(LIBS ${LIBS} ${PcapPlusPlus_LIBRARIES})
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


if(NOT WIN32)
#CoinUtils
if (CoinUtils)
message(STATUS "Enable CoinUtils")
add_definitions(-DUSE_CoinUtils)
find_package(CoinUtils QUIET)
if (CoinUtils_FOUND)
 message("-- Found CoinUtils: ${CoinUtils_INCLUDE_DIRS}")
 include_directories(${CoinUtils_INCLUDE_DIRS})
 set(LIBS ${LIBS} ${CoinUtils_LIBRARIES})
else(CoinUtils_FOUND)
 set(ADD_CoinUtils TRUE)
 message(STATUS "Will download CoinUtils library from https://github.com/coin-or/CoinUtils")
 include(ExternalCoinUtils)
 if(APPLE)
  set(CoinUtils_LIBRARY ${THIRDPARTY_INSTALL_PATH}/Install/CoinUtils/build/lib/libCoinUtils.dylib)
 elseif(UNIX)
  set(CoinUtils_LIBRARY ${THIRDPARTY_INSTALL_PATH}/Install/CoinUtils/build/lib/libCoinUtils.so)
 endif(APPLE)
endif()
endif()

#MP
if (MP)
message(STATUS "Enable AMPL MP")
add_definitions(-DUSE_MP)
find_package(MP QUIET)
if (MP_FOUND)
 message("-- Found AMPL MP: ${MP_INCLUDE_DIRS}")
 include_directories(${MP_INCLUDE_DIRS}/..)
 set(LIBS ${LIBS} ${MP_LIBRARIES})
else(MP_FOUND)
 set(ADD_MP TRUE)
 message(STATUS "Will download AMPL MP library from https://github.com/ampl/mp.git")
 include(ExternalMP)
 set(MP_LIBRARY ${THIRDPARTY_INSTALL_PATH}/Install/MP/build/lib/libmp.a)
endif()
endif()

# Find HiGHS?
if (HiGHS)
    message(STATUS "Enable HiGHS")
    add_definitions(-DUSE_HiGHS)
    find_package(HiGHS)
if (HiGHS_FOUND)
    message("-- Found HiGHS: ${HiGHS_INCLUDE_DIRS}")
    include_directories(${HiGHS_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${HiGHS_LIBRARIES})
else(HiGHS_FOUND)
    set(ADD_HiGHS TRUE)
    message(STATUS "Will download HiGHS library from https://github.com/ERGO-Code/HiGHS")
    include(ExternalHiGHS)
    if(APPLE)
        set(HiGHS_LIBRARIES ${THIRDPARTY_INSTALL_PATH}/Install/HiGHS/build/lib/libhighs.dylib)
        elseif(UNIX)
            set(HiGHS_LIBRARIES ${THIRDPARTY_INSTALL_PATH}/Install/HiGHS/build/lib/libhighs.so)
        elseif(WIN32)
            set(HiGHS_LIBRARIES ${THIRDPARTY_INSTALL_PATH}/Install/HiGHS/build/lib/libhighs.lib)
    endif(APPLE)
endif(HiGHS_FOUND)
endif(HiGHS)

endif()
if(Ipopt)
# Add Ipopt include dir and link to libraries
if(WIN32)
# include_directories("${PROJECT_SOURCE_DIR}/thirdparty/Ipopt")
include_directories("${PROJECT_SOURCE_DIR}/thirdparty/Ipopt/include/coin-or")
# add_subdirectory(thirdparty/Ipopt)
# find_library(IPOPT_LIBRARIES 
#     ipopt.dll
# 	HINTS /usr/local/lib
# 	HINTS ${PROJECT_SOURCE_DIR}/thirdparty/Ipopt
# 	HINTS ${IPOPT_ROOT_DIR}/lib
# )
endif(WIN32)

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
	set(IPOPT_LIBRARIES ${PROJECT_SOURCE_DIR}/thirdparty/Ipopt/libipopt.so)
elseif(WIN32) 
 set(IPOPT_LIBRARIES ${PROJECT_SOURCE_DIR}/thirdparty/Ipopt/libipopt-3.lib)
endif(APPLE)
include(ExternalIpopt)
endif (IPOPT_FOUND)
endif()

# Eigen
if(!WIN32)
find_package(Eigen QUIET)
if (EIGEN_FOUND)
 message("-- Found Eigen version ${EIGEN_VERSION}: ${EIGEN_INCLUDE_DIRS}")
else(EIGEN_FOUND)
 set(EIGEN_VERSION 3.3.7)
 include(ExternalEigen)
 set(EIGEN_INCLUDE_DIR ${THIRDPARTY_INSTALL_PATH}/Install/eigen)
 set(ADD_EIGEN TRUE)
endif (EIGEN_FOUND)
include_directories(${EIGEN_INCLUDE_DIR})
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

# Include Power dirs
include_directories(${PROJECT_SOURCE_DIR}/examples/Optimization/NonLinear/Power)
include_directories(${PROJECT_SOURCE_DIR}/examples/Optimization/Linear/Rostering)
# include_directories(${PROJECT_SOURCE_DIR}/solvers)

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

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-absolute-value")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath,'\$ORIGIN/lib'")


if(UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
  set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -pthread")
endif()

set(CMAKE_XCODE_ATTRIBUTE_LD_RUNPATH_SEARCH_PATHS "$ENV{DYLD_LIBRARY_PATH}")
add_subdirectory(src)
add_subdirectory(examples)
# add_subdirectory(solvers)

#Get user defined configuration including passing project source directory into c++ code
configure_file(GravityConfig.h.in ${PROJECT_SOURCE_DIR}/include/gravity/GravityConfig.h)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/gravity DESTINATION include)
