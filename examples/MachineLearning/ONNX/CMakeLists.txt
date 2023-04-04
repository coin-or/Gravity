find_package(Protobuf REQUIRED)

file(DOWNLOAD "https://raw.githubusercontent.com/onnx/onnx/master/onnx/onnx.proto3" "${CMAKE_CURRENT_BINARY_DIR}/onnx.proto")
file(DOWNLOAD "https://github.com/onnx/models/raw/main/vision/classification/mnist/model/mnist-12.onnx" "${PROJECT_BINARY_DIR}/mnist.onnx")
PROTOBUF_GENERATE_CPP(ONNX_PROTO_SRCS ONNX_PROTO_HDRS ${CMAKE_CURRENT_BINARY_DIR}/onnx.proto)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(
    onnx_parser
	ONNX/ONNX.cpp
    ${ONNX_PROTO_SRCS}
    ${ONNX_PROTO_HDRS}
)
message(${ONNX_PROTO_SRCS})
message(${ONNX_PROTO_HDRS})
target_link_libraries(onnx_parser ${PROTOBUF_LIBRARIES} gravity)