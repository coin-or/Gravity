#include(ExternalProject)

include_directories(${PROJECT_INCLUDE_DIR})

#message(STATUS "Adding Option Parser")
#ExternalProject_Add(opt_parser
#        SOURCE_DIR ${PROJECT_SOURCE_DIR}/third_party/cpp_option_parser
#        BUILD_IN_SOURCE 1
#        BUILD_COMMAND cmake .
#        INSTALL_COMMAND make
#)
#ExternalProject_Get_Property(opt_parser install_dir)
#add_library(option_parser STATIC IMPORTED)
#set_property(TARGET option_parser PROPERTY IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/third_party/cpp_option_parser/liboption_parser.a)
#add_dependencies(option_parser opt_parser)


# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
#SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
#SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:'$ORIGIN/lib'")

file(GLOB HEADERS ${PROJECT_INCLUDE_DIR} *.h)
file(GLOB SOURCES . *.cpp ${PROJECT_INCLUDE_DIR}/gravity/*.h)

if (NOT Ipopt)
  list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/IpoptProgram.cpp)
endif()
if (NOT Bonmin)
  list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/BonminProgram.cpp)
endif()
if (NOT Gurobi)
  list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/GurobiProgram.cpp)
endif()

if (NOT Cplex)
    list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/CplexProgram.cpp)
endif()

if (NOT Clp)
    list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/ClpProgram.cpp)
endif()

if (NOT Sdpa)
    list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/SdpaProgram.cpp)
endif()

if (NOT Mosek)
    list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/MosekProgram.cpp)
endif()


add_library(gravity STATIC ${SOURCES} ${HEADERS})
if(ADD_EIGEN)
 add_dependencies(gravity eigen)
endif(ADD_EIGEN)
add_dependencies(gravity opt_parser)
target_link_libraries(gravity ${CPLEX_LIBRARIES} ${GUROBI_LIBRARIES} ${IPOPT_LIBRARIES} ${BONMIN_LIBRARIES} ${SDPA_LIBRARIES} ${MOSEK_LIBRARIES} ${XLNT_LIBRARIES} ${OpenMPI_LIBRARIES} liboption_parser.a)
if(Cov)
  set_target_properties(gravity
  PROPERTIES
    COMPILE_FLAGS "-ftest-coverage -fprofile-arcs"
    LINK_FLAGS    "-lgcov --coverage -fprofile-arcs"
  )
endif()

SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
SET(CMAKE_INSTALL_RPATH "${PROJ_LIB_DIR}")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
