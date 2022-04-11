# Search supplied hint directories first if supplied.
set(TEASER_ROOT_DIR "$ENV{TEASER_ROOT_DIR}" CACHE PATH "TEASER root directory.")
if(TEASER_ROOT_DIR)
 message("Looking for Teaser in ${TEASER_ROOT_DIR}")
else(TEASER_ROOT_DIR)
 message("TEASER_ROOT_DIR not provided.")
endif(TEASER_ROOT_DIR)

find_path(TEASER_INCLUDE_DIR
  HINTS ${TEASER_ROOT_DIR}/teaser/include/teaser
)

find_library(TEASER_LIBRARY 
  libteaser_registration.so
  HINTS ${TEASER_ROOT_DIR}/build/teaser
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TEASER DEFAULT_MSG TEASER_LIBRARY TEASER_INCLUDE_DIR)



if(TEASER_FOUND)
	message("—- Found TEASER under ${TEASER_INCLUDE_DIR}")
    set(TEASER_INCLUDE_DIRS ${TEASER_INCLUDE_DIR})
    set(TEASER_LIBRARIES ${TEASER_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(TEASER_LIBRARIES "${TEASER_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(TEASER_FOUND)

mark_as_advanced(TEASER_LIBRARY TEASER_INCLUDE_DIR)