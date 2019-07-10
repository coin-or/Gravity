set(OpenMPI_ROOT_DIR "$ENV{OpenMPI_ROOT_DIR}" CACHE PATH "OpenMPI root directory.")
message("Looking for OpenMPI in ${OpenMPI_ROOT_DIR}")


find_path(OpenMPI_INCLUDE_DIR
	NAMES mpi.h 
	HINTS ${PROJECT_SOURCE_DIR}/third_party/OpenMPI/build/include
	HINTS /usr/local/include
	HINTS ${OpenMPI_ROOT_DIR}/include
	HINTS ${PROJECT_SOURCE_DIR}/include
)

if(APPLE)
find_library(OpenMPI_LIBRARY 
	libmpi.dylib
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/OpenMPI/build/lib
	HINTS ${OpenMPI_ROOT_DIR}/lib
	HINTS ${PROJECT_SOURCE_DIR}/lib
)
elseif(UNIX)
find_library(OpenMPI_LIBRARY 
	libmpi.so
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/OpenMPI/build/lib
	HINTS ${OpenMPI_ROOT_DIR}/lib
	HINTS ${PROJECT_SOURCE_DIR}/lib
)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenMPI DEFAULT_MSG OpenMPI_LIBRARY OpenMPI_INCLUDE_DIR)

if(OpenMPI_FOUND)
	message("—- Found OpenMPI under ${OpenMPI_INCLUDE_DIR}")
    set(OpenMPI_INCLUDE_DIRS ${OpenMPI_INCLUDE_DIR})
    set(OpenMPI_LIBRARIES ${OpenMPI_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(OpenMPI_LIBRARIES "${OpenMPI_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(OpenMPI_FOUND)

mark_as_advanced(OpenMPI_LIBRARY OpenMPI_INCLUDE_DIR)
