set(ARMADILLO_ROOT_DIR "$ENV{ARMADILLO_ROOT_DIR}" CACHE PATH "ARMADILLO root directory.")

find_path(ARMADILLO_INCLUDE_DIR
        NAMES armadillo
        HINTS ARMADILLO_ROOT_DIR/usr/include/
        HINTS /usr/include
        HINTS /usr/local/include/
        )

find_library(ARMADILLO_LIBRARY
        libarmadillo.so
        HINTS ARMADILLO_ROOT_DIR/usr/lib/x86_64-linux-gnu/
        HINTS /usr/lib/x86_64-linux-gnu/
        )

find_library(LAPACK_LIBRARY
        liblapack.so
        HINTS /usr/lib/
        )

find_library(BLAS_LIBRARY
        libblas.so
        HINTS /usr/lib/
        )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARMADILLO DEFAULT_MSG ARMADILLO_LIBRARY ARMADILLO_INCLUDE_DIR LAPACK_LIBRARY BLAS_LIBRARY)

if(ARMADILLO_FOUND)
    message("—- Found Armadillo under ${ARMADILLO_INCLUDE_DIR}")
    message("—- Library ${ARMADILLO_LIBRARY} ${LAPACK_LIBRARY} ${BLAS_LIBRARY}")
    set(ARMADILLO_INCLUDE_DIRS ${ARMADILLO_INCLUDE_DIR})
    set(ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARY})
endif(ARMADILLO_FOUND)