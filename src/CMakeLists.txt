include_directories(${PROJECT_INCLUDE_DIR})

file(GLOB HEADERS ${PROJECT_INCLUDE_DIR} *.h)
file(GLOB SOURCES . *.cpp ${PROJECT_INCLUDE_DIR}/gravity/*.h)

add_library(gravity ${SOURCES} ${HEADERS})
target_link_libraries(gravity ${GUROBI_LIBRARIES} ${IPOPT_LIBRARIES} ${BONMIN_LIBRARIES})

add_executable(gravity_test main.cpp)
target_link_libraries(gravity_test gravity)
