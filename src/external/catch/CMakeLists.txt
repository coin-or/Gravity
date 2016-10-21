cmake_minimum_required(VERSION 2.8.8)
project(catch_builder CXX)
include(ExternalProject)
include_directories(${PROJECT_INCLUDE_DIR})
find_package(Git REQUIRED)

ExternalProject_Add(
    catch
    PREFIX ${CMAKE_BINARY_DIR}/catch
    GIT_REPOSITORY https://github.com/philsquared/Catch.git
    TIMEOUT 10
    UPDATE_COMMAND ${GIT_EXECUTABLE} pull
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    LOG_DOWNLOAD ON
   )

# Expose required variable (CATCH_INCLUDE_DIR) to parent scope
ExternalProject_Get_Property(catch source_dir)
set(CATCH_INCLUDE_DIR ${source_dir}/include CACHE INTERNAL "Path to include folder for Catch")

add_executable(UnitTest catch.hpp unitTest.cpp catch.cpp)
target_link_libraries(UnitTest gravity)
