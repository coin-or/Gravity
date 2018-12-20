include_directories(${PROJECT_INCLUDE_DIR})





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

set(CPP_OPTION_PARSER option_parser)

add_library(gravity STATIC ${SOURCES} ${HEADERS})
target_link_libraries(gravity ${CPLEX_LIBRARIES} ${GUROBI_LIBRARIES} ${IPOPT_LIBRARIES} ${BONMIN_LIBRARIES} ${SDPA_LIBRARIES} ${MOSEK_LIBRARIES} ${XLNT_LIBRARIES} ${CPP_OPTION_PARSER})
set_target_properties(gravity
PROPERTIES
    COMPILE_FLAGS "-ftest-coverage -fprofile-arcs"
    LINK_FLAGS    "-lgcov --coverage"
)
