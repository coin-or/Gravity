set(H5CPP_ROOT_DIR "$ENV{H5CPP_ROOT_DIR}" CACHE PATH "H5CPP root directory.")
message("Looking for H5CPP in ${H5CPP_ROOT_DIR}")


include(FindPackageHandleStandardArgs)

find_path(H5CPP_INCLUDE_DIR
	NAMES hdf5.hpp 
	HINTS /usr/local/include
	HINTS /usr/local/include/h5cpp
	HINTS ${H5CPP_ROOT_DIR}/include
	HINTS ${H5CPP_ROOT_DIR}/src/h5cpp
	HINTS ${PROJECT_SOURCE_DIR}/thirdparty/H5CPP/src
)


find_path(HDF5_INCLUDE_DIR
        NAMES hdf5.h
	HINTS /usr/local/opt/hdf5/include
        HINTS ${HDF5_ROOT_DIR}/src/
        HINTS ${PROJECT_SOURCE_DIR}/thirdparty/HDF5/src
        HINTS /usr/local/include
)
include_directories(${HDF5_INCLUDE_DIR})
if(WIN32)


find_library(H5CPP_LIBRARY 
	libh5cpp.lib
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/thirdparty/H5CPP/lib
	HINTS ${H5CPP_ROOT_DIR}/lib
	HINTS ${H5CPP_ROOT_DIR}/build/lib
)

find_library(HDF5_LIBRARY
        libhdf5.lib
        HINTS /usr/local/lib
        HINTS ${PROJECT_SOURCE_DIR}/thirdparty/H5CPP/lib
        HINTS ${H5CPP_ROOT_DIR}/lib
        HINTS ${H5CPP_ROOT_DIR}/build/lib
)

find_package_handle_standard_args(H5CPP DEFAULT_MSG H5CPP_LIBRARY H5CPP_INCLUDE_DIR HDF5_INCLUDE_DIR)
elseif(APPLE)

find_library(H5CPP_LIBRARY 
	libh5cpp.dylib	
	HINTS ${H5CPP_ROOT_DIR}/lib
	HINTS ${H5CPP_ROOT_DIR}/build/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/H5CPP/build/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/H5CPP/lib	
	HINTS /usr/local/lib
)

find_library(HDF5_LIBRARY
        libhdf5.dylib
	HINTS /usr/local/opt/hdf5/lib
        HINTS ${PROJECT_SOURCE_DIR}/thirdparty/H5CPP/lib
        HINTS ${H5CPP_ROOT_DIR}/lib
        HINTS ${H5CPP_ROOT_DIR}/build/lib
        HINTS /usr/local/lib
)

find_package_handle_standard_args(H5CPP DEFAULT_MSG H5CPP_LIBRARY HDF5_LIBRARY H5CPP_INCLUDE_DIR HDF5_INCLUDE_DIR)
elseif(UNIX)

find_library(H5CPP_LIBRARY 
	libh5cpp.so	
	HINTS ${H5CPP_ROOT_DIR}/lib
	HINTS ${H5CPP_ROOT_DIR}/build/lib
	HINTS ${PROJECT_SOURCE_DIR}/thirdparty/H5CPP/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/H5CPP/build/lib
	HINTS /usr/local/lib
)

find_library(HDF5_LIBRARY
        libhdf5.so
        HINTS ${H5CPP_ROOT_DIR}/lib
        HINTS ${H5CPP_ROOT_DIR}/build/lib
        HINTS ${PROJECT_SOURCE_DIR}/thirdparty/H5CPP/lib
        HINTS ${PROJECT_SOURCE_DIR}/third_party/H5CPP/build/lib
        HINTS /usr/local/lib
)

find_package_handle_standard_args(H5CPP DEFAULT_MSG H5CPP_LIBRARY HDF5_LIBRARY H5CPP_INCLUDE_DIR HDF5_INCLUDE_DIR)
endif(WIN32)

if(H5CPP_FOUND)
	message("—- Found H5CPP include dir under ${H5CPP_INCLUDE_DIR}")
	message("—- Found HDF5 include dir under ${HDF5_INCLUDE_DIR}")
	message("—- Found H5CPP lib at ${H5CPP_LIBRARY}")
	message("—- Found HDF5 lib at ${HDF5_LIBRARY}")
    set(H5CPP_INCLUDE_DIRS ${H5CPP_INCLUDE_DIR})
#    set(H5CPP_INCLUDE_DIRS ${H5CPP_INCLUDE_DIR} ${HDF5_INCLUDE_DIR})
    set(H5CPP_LIBRARIES ${H5CPP_LIBRARY} ${HDF5_LIBRARY})
else (H5CPP_FOUND)
 message("Cannot find H5CPP, will try pulling it from github.")
endif(H5CPP_FOUND)

mark_as_advanced(H5CPP_LIBRARY H5CPP_INCLUDE_DIR)

