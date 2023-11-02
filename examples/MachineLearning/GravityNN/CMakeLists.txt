include(${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/cmake/CPM.cmake)

include(FetchContent)
CPMAddPackage("gh:CLIUtils/CLI11#v2.3.2")
CPMAddPackage("gh:onnx/onnx#v1.15.0")

file(GLOB NN_SOURCES 
  ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network/*.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network/Layers/*.hpp
)

add_executable(
  GravityNN
	MachineLearning/GravityNN/GravityNN.cpp
	${ONNX_PROTO_SRCS}	
	${ONNX_PROTO_HDRS}
	${NN_SOURCES}
)

target_include_directories(GravityNN PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network  ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network/Layers ${ONNX_PROTO_HDRS} ${cli11_SOURCE_DIR}/include)
target_link_libraries(GravityNN gravity CLI11 onnx)

find_package(OpenMP)
add_executable(
  tighten
	MachineLearning/GravityNN/tighten.cpp
	${ONNX_PROTO_SRCS}
	${ONNX_PROTO_HDRS}
	${NN_SOURCES}
)

target_include_directories(tighten PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network  ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network/Layers ${ONNX_PROTO_HDRS} ${cli11_SOURCE_DIR}/include)
target_link_libraries(tighten gravity CLI11 onnx)
if (OpenMPI)
	include_directories(${MPI_INCLUDE_PATH})
	target_link_libraries(tighten ${MPI_LIBRARIES})
endif(OpenMPI)
