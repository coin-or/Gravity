set(HDF5_ROOT_DIR "$ENV{HDF5_ROOT_DIR}" CACHE PATH "HDF5 root directory.")
message("Looking for HDF5 in ${HDF5_ROOT_DIR}")


include(FindPackageHandleStandardArgs)


find_path(HDF5_INCLUDE_DIR
        NAMES hdf5.h
	HINTS /usr/local/opt/hdf5/include
        HINTS ${HDF5_ROOT_DIR}/src/
        HINTS ${PROJECT_SOURCE_DIR}/thirdparty/HDF5/src
        HINTS /usr/local/include
)
include_directories(${HDF5_INCLUDE_DIR})
if(WIN32)


find_library(HDF5_LIBRARY
        libhdf5.lib
        HINTS /usr/local/lib
        HINTS ${PROJECT_SOURCE_DIR}/thirdparty/HDF5/lib
        HINTS ${HDF5_ROOT_DIR}/lib
        HINTS ${HDF5_ROOT_DIR}/build/lib
)

elseif(APPLE)


find_library(HDF5_LIBRARY
        libhdf5.dylib
	HINTS /usr/local/opt/hdf5/lib
        HINTS ${PROJECT_SOURCE_DIR}/thirdparty/HDF5/lib
        HINTS ${HDF5_ROOT_DIR}/lib
        HINTS ${HDF5_ROOT_DIR}/build/lib
        HINTS /usr/local/lib
)

elseif(UNIX)

find_library(HDF5_LIBRARY
        libhdf5.so
        HINTS ${HDF5_ROOT_DIR}/lib
        HINTS ${HDF5_ROOT_DIR}/build/lib
        HINTS ${PROJECT_SOURCE_DIR}/thirdparty/HDF5/lib
        HINTS ${PROJECT_SOURCE_DIR}/third_party/HDF5/build/lib
        HINTS /usr/local/lib
)


endif(WIN32)

find_package_handle_standard_args(HDF5 DEFAULT_MSG HDF5_LIBRARY HDF5_INCLUDE_DIR)

if(HDF5_FOUND)
	message("—- Found HDF5 include dir under ${HDF5_INCLUDE_DIR}")
	message("—- Found HDF5 lib at ${HDF5_LIBRARY}")
    set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR})
#    set(HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIR} ${HDF5_INCLUDE_DIR})
    set(HDF5_LIBRARIES ${HDF5_LIBRARY})
else (HDF5_FOUND)
 message("Cannot find HDF5, please install it (e.g., brew install hdf5).")
endif(HDF5_FOUND)

mark_as_advanced(HDF5_LIBRARY HDF5_INCLUDE_DIR)

