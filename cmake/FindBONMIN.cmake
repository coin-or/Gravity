set(BONMIN_ROOT_DIR "/home/kbestuzheva/Bonmin-1.8.4" CACHE PATH "BONMIN root directory.")

message("Looking for Bonmin in ${BONMIN_ROOT_DIR}")

#string(REGEX MATCH "[0-9]+" BONMIN_VERSION "${BONMIN_ROOT_DIR}")

#message("Ipopt version ${BONMIN_VERSION}")

#string(SUBSTRING ${BONMIN_VERSION} 0 2 BONMIN_VERSION_SHORT)


find_path(BONMIN_INCLUDE_DIR
        NAMES BonTMINLP.hpp
        HINTS /usr/local/include/coin
        HINTS ${BONMIN_ROOT_DIR}/include/coin
        HINTS ${BONMIN_ROOT_DIR}/include
)

if(APPLE)
    find_library(BONMIN_LIBRARY
            libbonmin.dylib
            HINTS /usr/local/lib
            HINTS ${BONMIN_ROOT_DIR}/lib
    )
elseif(UNIX)
    find_library(BONMIN_LIBRARY
            libbonmin.so
            HINTS /usr/local/lib
            HINTS ${BONMIN_ROOT_DIR}/lib
    )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BONMIN DEFAULT_MSG BONMIN_LIBRARY BONMIN_INCLUDE_DIR)

if(BONMIN_FOUND)
    message("—- Found Bonmin under ${BONMIN_INCLUDE_DIR}")
    set(BONMIN_INCLUDE_DIRS ${BONMIN_INCLUDE_DIR})
    set(BONMIN_LIBRARIES ${BONMIN_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(BONMIN_LIBRARIES "${BONMIN_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(BONMIN_FOUND)

mark_as_advanced(BONMIN_LIBRARY BONMIN_INCLUDE_DIR)
