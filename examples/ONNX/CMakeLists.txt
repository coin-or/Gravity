#############################
# CPM
#############################
set(CPM_DOWNLOAD_VERSION 0.38.0)

if(CPM_SOURCE_CACHE)
  set(CPM_DOWNLOAD_LOCATION "${CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
elseif(DEFINED ENV{CPM_SOURCE_CACHE})
  set(CPM_DOWNLOAD_LOCATION "$ENV{CPM_SOURCE_CACHE}/cpm/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
else()
  set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_DOWNLOAD_VERSION}.cmake")
endif()

# Expand relative path. This is important if the provided path contains a tilde (~)
get_filename_component(CPM_DOWNLOAD_LOCATION ${CPM_DOWNLOAD_LOCATION} ABSOLUTE)

function(download_cpm)
  message(STATUS "Downloading CPM.cmake to ${CPM_DOWNLOAD_LOCATION}")
  file(DOWNLOAD
       https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_DOWNLOAD_VERSION}/CPM.cmake
       ${CPM_DOWNLOAD_LOCATION}
  )
endfunction()

if(NOT (EXISTS ${CPM_DOWNLOAD_LOCATION}))
  download_cpm()
else()
  # resume download if it previously failed
  file(READ ${CPM_DOWNLOAD_LOCATION} check)
  if("${check}" STREQUAL "")
    download_cpm()
  endif()
endif()

include(${CPM_DOWNLOAD_LOCATION})

####################################
# ONNX
####################################
CPMAddPackage("gh:gabime/spdlog#v1.11.0")
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
target_link_libraries(onnx_parser ${PROTOBUF_LIBRARIES} gravity spdlog::spdlog)