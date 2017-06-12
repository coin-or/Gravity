include_directories(${PROJECT_INCLUDE_DIR})

add_executable(gravity_test ${PROJECT_SOURCE_DIR}/examples/gravity101.cpp)
add_executable(stable_set ${PROJECT_SOURCE_DIR}/examples/Stable_Set.cpp)
target_link_libraries(stable_set gravity)
target_link_libraries(gravity_test gravity)
