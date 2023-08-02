find_package(Protobuf REQUIRED)

file(DOWNLOAD "https://raw.githubusercontent.com/onnx/onnx/master/onnx/onnx.proto3" "${CMAKE_CURRENT_BINARY_DIR}/onnx.proto")
PROTOBUF_GENERATE_CPP(ONNX_PROTO_SRCS ONNX_PROTO_HDRS ${CMAKE_CURRENT_BINARY_DIR}/onnx.proto)
add_custom_target(proto_dep DEPENDS ${ONNX_PROTO_SRCS} ${ONNX_PROTO_HDRS})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

find_package(PkgConfig QUIET REQUIRED)
pkg_check_modules(PC_Protobuf protobuf)
message(${PC_Protobuf_INCLUDEDIR})

include_directories(${PC_Protobuf_INCLUDEDIR} ${ONNX_PROTO_HDRS} ${ONNX_PROTO_SRCS})

include(FetchContent)
FetchContent_Declare(
  cli11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11
  GIT_TAG        v2.3.2
)

file(GLOB NN_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network/*.hpp ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network/Layers/*.hpp)
FetchContent_MakeAvailable(cli11)
message(${cli11_SOURCE_DIR})
add_executable(
  GravityNN
	MachineLearning/GravityNN/GravityNN.cpp
	${ONNX_PROTO_SRCS}	
	${ONNX_PROTO_HDRS}
	${NN_SOURCES}
)

add_dependencies(GravityNN proto_dep)


target_include_directories(GravityNN PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network  ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network/Layers ${ONNX_PROTO_HDRS} ${cli11_SOURCE_DIR}/include)
target_link_libraries(GravityNN ${PROTOBUF_LIBRARIES} gravity CLI11)

find_package(OpenMP)
add_executable(
  tighten
	MachineLearning/GravityNN/tighten.cpp
        ${ONNX_PROTO_SRCS}
        ${ONNX_PROTO_HDRS}
        ${NN_SOURCES}
)
add_dependencies(tighten proto_dep)
target_include_directories(tighten PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network  ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network/Layers ${ONNX_PROTO_HDRS} ${cli11_SOURCE_DIR}/include)
target_link_libraries(tighten ${PROTOBUF_LIBRARIES} gravity CLI11)
