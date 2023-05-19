find_package(Protobuf REQUIRED)

file(DOWNLOAD "https://raw.githubusercontent.com/onnx/onnx/master/onnx/onnx.proto3" "${CMAKE_CURRENT_BINARY_DIR}/onnx.proto")
PROTOBUF_GENERATE_CPP(ONNX_PROTO_SRCS ONNX_PROTO_HDRS ${CMAKE_CURRENT_BINARY_DIR}/onnx.proto)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

find_package(PkgConfig QUIET REQUIRED)
pkg_check_modules(PC_Protobuf protobuf)
message(${PC_Protobuf_INCLUDEDIR})

include_directories(${PC_Protobuf_INCLUDEDIR})


file(GLOB GravityNN_HEADERS MachineLearning/GravityNN/*.hpp)
add_executable(
  GravityNN
	MachineLearning/GravityNN/GravityNN.cpp
	${GravityNN_HEADERS}
	${ONNX_PROTO_SRCS}
  ${ONNX_PROTO_HDRS}
)
target_link_libraries(GravityNN ${PROTOBUF_LIBRARIES} gravity)

add_executable(
  tighten
	MachineLearning/GravityNN/tighten.cpp
	${GravityNN_HEADERS}
	${ONNX_PROTO_SRCS}
  ${ONNX_PROTO_HDRS}
)
target_link_libraries(tighten ${PROTOBUF_LIBRARIES} gravity)