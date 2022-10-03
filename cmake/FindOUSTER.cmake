set(OUSTER_ROOT_DIR "$ENV{OUSTER_ROOT_DIR}" CACHE PATH "OUSTER root directory.")
message("Looking for OUSTER in ${OUSTER_ROOT_DIR}")


find_path(OUSTER_INCLUDE_DIR
    NAMES lidar_scan.h 
    HINTS /usr/local/include/ouster
    HINTS ${OUSTER_ROOT_DIR}/include/ouster
    HINTS ${PROJECT_SOURCE_DIR}/third_party/OUSTER/build/include/coin
)

find_path(OUSTER_INCLUDE_DIR2
    NAMES optional.hpp
    HINTS /usr/local/include/optional-lite/nonstd/
    HINTS ${OUSTER_ROOT_DIR}/include/ouster
    HINTS ${PROJECT_SOURCE_DIR}/third_party/OUSTER/build/include/coin
)

find_library(OUSTER_LIBRARY 
    libouster_pcap.a
    HINTS /usr/local/lib
    HINTS ${PROJECT_SOURCE_DIR}/third_party/OUSTER/build/lib
    HINTS ${OUSTER_ROOT_DIR}/build/lib
)

find_library(OUSTER_LIBRARY2 
    libouster_client.a
    HINTS /usr/local/lib
    HINTS ${PROJECT_SOURCE_DIR}/third_party/OUSTER/build/lib
    HINTS ${OUSTER_ROOT_DIR}/build/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OUSTER DEFAULT_MSG OUSTER_LIBRARY OUSTER_INCLUDE_DIR)

if(OUSTER_FOUND)
    message("—- Found OUSTER under ${OUSTER_INCLUDE_DIR} and ${OUSTER_INCLUDE_DIR2}")
    include_directories(${OUSTER_INCLUDE_DIR}/../)
    include_directories(${OUSTER_INCLUDE_DIR2}/../)
    set(OUSTER_INCLUDE_DIRS ${OUSTER_INCLUDE_DIR} ${OUSTER_INCLUDE_DIR2})
    set(OUSTER_LIBRARIES ${OUSTER_LIBRARY} ${OUSTER_LIBRARY2})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(OUSTER_LIBRARIES "${OUSTER_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
else (OUSTER_FOUND)
 message("Cannot find OUSTER, will try pulling it from github.")
endif(OUSTER_FOUND)

mark_as_advanced(OUSTER_LIBRARY OUSTER_INCLUDE_DIR)
