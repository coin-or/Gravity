include_directories(${PROJECT_INCLUDE_DIR})

file(GLOB HEADERS ${PROJECT_INCLUDE_DIR} *.h)
file(GLOB SOURCES . *.cpp ${PROJECT_INCLUDE_DIR}/gravity/*.h)

if (NOT ENABLE_IPOPT)
  list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/IpoptProgram.cpp)
endif()
if (NOT ENABLE_BONMIN)
  list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/BonminProgram.cpp)
endif()
if (NOT ENABLE_GUROBI)
  list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/GurobiProgram.cpp)
endif()

if (NOT ENABLE_CPLEX)
    list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/CplexProgram.cpp)
endif()

if (NOT ENABLE_SDPA)
    list(REMOVE_ITEM SOURCES ${PROJECT_SOURCE_DIR}/src/SdpaProgram.cpp)
endif()

add_library(gravity ${SOURCES} ${HEADERS})
target_link_libraries(gravity ${CPLEX_LIBRARIES} ${GUROBI_LIBRARIES}  ${IPOPT_LIBRARIES} ${BONMIN_LIBRARIES} ${SDPA_LIBRARIES})

