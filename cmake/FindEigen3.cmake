set(EIGEN3_ROOT_DIR "$ENV{EIGEN3_ROOT_DIR}" CACHE PATH "EIGEN3 root directory.")
message("Looking for EIGEN3 in ${EIGEN3_ROOT_DIR}")

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DEIGEN3_ROOT_DIR:PATH=${EIGEN3_ROOT_DIR}")
include_directories(${EIGEN3_INCLUDE_DIRS})

set(EIGEN3_INCLUDE_DIRS ${EIGEN3_ROOT_DIR})

