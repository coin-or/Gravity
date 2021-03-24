include_directories(${PROJECT_INCLUDE_DIR})

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
#SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


# the RPATH to be used when installing, but only if it's not a system directory
#SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib:'$ORIGIN/lib'")
file(GLOB ARMO_HEADERS ${PROJECT_SOURCE_DIR}/ARMO/*.hpp)
file(GLOB ARMO_SOURCES ${PROJECT_SOURCE_DIR}/ARMO/*.cpp)


if (NOT Matplot)
  list(REMOVE_ITEM HEADERS ${PROJECT_INCLUDE_DIR}/matplotlibcpp.h)
endif()

add_executable(ARMO ${ARMO_SOURCES} ${ARMO_HEADERS})
target_link_libraries(ARMO gravity)
