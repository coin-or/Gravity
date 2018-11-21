file(GLOB dirs ${PROJECT_SOURCE_DIR}/third_party/armadillo-*) 
foreach(d in ${dirs})
   string(REGEX MATCH "[0-9]+.[0-9]+.[0-9]+" ARMADILLO_VERSION "${d}")
   set(ARMADILLO_ROOT_DIR ${d})
endforeach(d)

message("Armadillo root dir ${ARMADILLO_ROOT_DIR}")
#string(SUBSTRING ${ARMADILLO_ROOT_DIR} 10 -1 ARMADILLO_VERSION)

message("Armadillo version ${ARMADILLO_VERSION}")


#find_path(ARMADILLO_INCLUDE_DIR
#        NAMES armadillo
#        HINTS ARMADILLO_ROOT_DIR/usr/include/
#        HINTS /usr/include
#        HINTS /usr/local/include/
#        )

#find_library(ARMADILLO_LIBRARY
#        libarmadillo.so
#        HINTS ARMADILLO_ROOT_DIR/usr/lib/x86_64-linux-gnu/
#        HINTS /usr/lib/x86_64-linux-gnu/
#        )

      #find_library(LAPACK_LIBRARY
      #  liblapack.so
      #  HINTS /usr/lib/
      #  )

      #find_library(BLAS_LIBRARY
      #  libblas.so
      #  HINTS /usr/lib/
      #  )

      #include(FindPackageHandleStandardArgs)
#find_package_handle_standard_args(ARMADILLO DEFAULT_MSG ARMADILLO_LIBRARY ARMADILLO_INCLUDE_DIR LAPACK_LIBRARY BLAS_LIBRARY)
#find_package_handle_standard_args(ARMADILLO DEFAULT_MSG ARMADILLO_LIBRARY ARMADILLO_INCLUDE_DIR)

#if(ARMADILLO_FOUND)
#    message("—- Found Armadillo under ${ARMADILLO_INCLUDE_DIR}")
#    message("—- Library ${ARMADILLO_LIBRARY}")
#    set(ARMADILLO_INCLUDE_DIRS ${ARMADILLO_INCLUDE_DIR})
#    set(ARMADILLO_LIBRARIES ${ARMADILLO_LIBRARY})
#endif(ARMADILLO_FOUND)
