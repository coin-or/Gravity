set(LASlib_ROOT_DIR "$ENV{LASlib_ROOT_DIR}" CACHE PATH "LASlib root directory.")
message("Looking for LASlib in ${LASlib_ROOT_DIR}")


find_path(LASlib_INCLUDE_DIR
	NAMES lasutility.hpp
	HINTS /usr/local/include/LASlib
	HINTS ${LASlib_ROOT_DIR}/include
	HINTS ${LASlib_ROOT_DIR}/inc
	HINTS ${LASlib_ROOT_DIR}
)

find_path(LASlib_INCLUDE_DIR2
        NAMES mydefs.hpp
        HINTS /usr/local/include/LASzip
        HINTS ${LASlib_ROOT_DIR}../LASzip/src
        HINTS ${LASlib_ROOT_DIR}/include
        HINTS ${LASlib_ROOT_DIR}/inc
)

find_library(LASlib_LIBRARY 
	liblas.a
	HINTS /usr/local/lib/LASlib
	HINTS ${LASlib_ROOT_DIR}/lib
	HINTS ${LASlib_ROOT_DIR}/bin
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LASlib DEFAULT_MSG LASlib_LIBRARY LASlib_INCLUDE_DIR LASlib_INCLUDE_DIR2)

if(LASlib_FOUND)
	message("—- Found LASlib under ${LASlib_INCLUDE_DIR}")
    set(LASlib_INCLUDE_DIRS ${LASlib_INCLUDE_DIR} ${LASlib_INCLUDE_DIR2})
    set(LASlib_LIBRARIES ${LASlib_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(LASlib_LIBRARIES "${LASlib_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(LASlib_FOUND)

mark_as_advanced(LASlib_LIBRARY LASlib_INCLUDE_DIR)
