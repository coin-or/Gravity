set(SMTLIB_ROOT_DIR "$ENV{SMTLIB_ROOT_DIR}" CACHE PATH "SMTLIB root directory.")
message("Looking for SMTLIB in ${SMTLIB_ROOT_DIR}")


find_path(SMTLIB_INCLUDE_DIR
	NAMES smtlib2utils.h 
	HINTS ${PROJECT_SOURCE_DIR}/third_party/SMTLIB/
	HINTS /usr/local/include	
	HINTS ${SMTLIB_ROOT_DIR}
)

find_library(SMTLIB_LIBRARY 
	libsmtlib2parser.a
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/SMTLIB/
	HINTS ${SMTLIB_ROOT_DIR}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SMTLIB DEFAULT_MSG SMTLIB_LIBRARY SMTLIB_INCLUDE_DIR)

if(SMTLIB_FOUND)
	message("—- Found SMTLIB under ${SMTLIB_INCLUDE_DIR}")
    set(SMTLIB_INCLUDE_DIRS ${SMTLIB_INCLUDE_DIR})
    set(SMTLIB_LIBRARIES ${SMTLIB_LIBRARY})
endif(SMTLIB_FOUND)

mark_as_advanced(SMTLIB_LIBRARY SMTLIB_INCLUDE_DIR)
