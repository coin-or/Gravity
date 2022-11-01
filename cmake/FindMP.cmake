set(MP_ROOT_DIR "$ENV{MP_ROOT_DIR}" CACHE PATH "MP root directory.")
message("Looking for MP in ${MP_ROOT_DIR}")


find_path(MP_INCLUDE_DIR
	NAMES nl.h 
	HINTS /usr/local/include/
	HINTS ${MP_ROOT_DIR}/include/mp
	HINTS ${PROJECT_SOURCE_DIR}/third_party/MP/build/include/coin
)

find_library(MP_LIBRARY 
	libmp.a
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/MP/build/lib
	HINTS ${MP_ROOT_DIR}/build/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MP DEFAULT_MSG MP_LIBRARY MP_INCLUDE_DIR)

if(MP_FOUND)
	message("—- Found MP under ${MP_INCLUDE_DIR}")
    set(MP_INCLUDE_DIRS ${MP_INCLUDE_DIR})
    set(MP_LIBRARIES ${MP_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(MP_LIBRARIES "${MP_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
else (MP_FOUND)
 message("Cannot find MP, will try pulling it from github.")
endif(MP_FOUND)

mark_as_advanced(MP_LIBRARY MP_INCLUDE_DIR)
