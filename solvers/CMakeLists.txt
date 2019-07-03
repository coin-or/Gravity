include_directories(${PROJECT_INCLUDE_DIR})



# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
#SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
#SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:'$ORIGIN/lib'")
file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/solvers/*.h)
file(GLOB HEADERS1 ${PROJECT_SOURCE_DIR}/examples/Optimization/NonLinear/Power/*.h)

file(GLOB ALPINEX_SOURCES ${PROJECT_SOURCE_DIR}/solvers/AlpineX/*.cpp)
file(GLOB Power_SOURCES ${PROJECT_SOURCE_DIR}/examples/Optimization/NonLinear/Power/*.cpp)

add_executable(AlpineX ${ALPINEX_SOURCES} ${Power_SOURCES} ${HEADERS} ${HEADERS1})
target_link_libraries(AlpineX gravity)
