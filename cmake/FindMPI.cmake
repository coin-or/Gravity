set(MPI_ROOT_DIR "$ENV{MPI_ROOT_DIR}" CACHE PATH "MPI root directory.")
message("Looking for MPI in ${MPI_ROOT_DIR}")


find_path(MPI_INCLUDE_DIR
	NAMES mpi.h 
	HINTS ${PROJECT_SOURCE_DIR}/third_party/MPI/build/include/coin
	HINTS /usr/local/include	
	HINTS ${MPI_ROOT_DIR}/include
)

if(APPLE)
find_library(MPI_LIBRARY 
	libmpi.dylib
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/MPI/build/lib
	HINTS ${MPI_ROOT_DIR}/lib
)
elseif(UNIX)
find_library(MPI_LIBRARY 
	libmpi.so
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/MPI/build/lib
	HINTS ${MPI_ROOT_DIR}/lib
)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPI DEFAULT_MSG MPI_LIBRARY MPI_INCLUDE_DIR)

if(MPI_FOUND)
	message("—- Found MPI under ${MPI_INCLUDE_DIR}")
    set(MPI_INCLUDE_DIRS ${MPI_INCLUDE_DIR})
    set(MPI_LIBRARIES ${MPI_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(MPI_LIBRARIES "${MPI_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(MPI_FOUND)

mark_as_advanced(MPI_LIBRARY MPI_INCLUDE_DIR)
