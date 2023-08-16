
find_package(PkgConfig QUIET REQUIRED)
pkg_check_modules(PC_Protobuf Protobuf)
find_package(Protobuf REQUIRED HINTS /usr/local/)
message(${PC_Protobuf_INCLUDEDIR})
message(${PC_Protobuf_LIBRARIES})


#file(DOWNLOAD "https://raw.githubusercontent.com/onnx/onnx/master/onnx/onnx.proto3" "${CMAKE_CURRENT_BINARY_DIR}/onnx.proto")

#protobuf_generate(
#    TARGET gravity
#    PROTOS ${CMAKE_CURRENT_BINARY_DIR}/onnx.proto
#    OUT_VAR PROTO_GENERATED_FILES
#    IMPORT_DIRS ${PROTO_IMPORT_DIRS}
#    PROTOC_OUT_DIR "${PROTO_BINARY_DIR}")
#set_source_files_properties(${PROTO_GENERATED_FILES} PROPERTIES SKIP_UNITY_BUILD_INCLUSION on)


#PROTOBUF_GENERATE_CPP(ONNX_PROTO_SRCS ONNX_PROTO_HDRS ${CMAKE_CURRENT_BINARY_DIR}/onnx.proto)
#add_custom_target(proto_dep DEPENDS ${ONNX_PROTO_SRCS} ${ONNX_PROTO_HDRS})
#include_directories(${CMAKE_CURRENT_BINARY_DIR})


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

find_library(ONNX_LIB
            libonnx.a
            HINTS /usr/local/lib
            HINTS $ENV{HOME}/Utils/onnx/build
    )

find_library(ONNX_PROTO_LIB
            libonnx_proto.a
            HINTS /usr/local/lib
            HINTS $ENV{HOME}/Utils/onnx/build
    )

find_library(PROTOBUF_LIB
        libprotobuf.a
        HINTS /usr/local/lib
        HINTS $ENV{HOME}/Utils/onnx/build
)

#add_dependencies(GravityNN proto_dep)


target_include_directories(GravityNN PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network  ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network/Layers ${ONNX_PROTO_HDRS} ${cli11_SOURCE_DIR}/include)
target_link_libraries(GravityNN ${PROTOBUF_LIB} gravity CLI11 ${ONNX_LIB} ${ONNX_PROTO_LIB})

find_package(OpenMP)
add_executable(
  tighten
	MachineLearning/GravityNN/tighten.cpp
        ${ONNX_PROTO_SRCS}
        ${ONNX_PROTO_HDRS}
        ${NN_SOURCES}
)
#target_link_libraries(tighten absl::base absl::synchronization absl::strings)
#target_link_libraries(tighten ${ABSL_BASE_LIBRARY})

#add_dependencies(tighten proto_dep)
target_include_directories(tighten PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network  ${CMAKE_CURRENT_SOURCE_DIR}/MachineLearning/GravityNN/network/Layers ${ONNX_PROTO_HDRS} ${cli11_SOURCE_DIR}/include)
target_link_libraries(tighten ${PROTOBUF_LIB} gravity CLI11 ${ONNX_LIB}  ${ONNX_PROTO_LIB})
#target_link_libraries(tighten mpi)
