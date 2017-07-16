set(MOSEK_ROOT_DIR "$ENV{MOSEK_ROOT_DIR}" CACHE PATH "MOSEK root directory.")
message("Looking for MOSEK in ${MOSEK_ROOT_DIR}")


find_path(MOSEK_INCLUDE_DIR	 fusion.h	HINTS "${MOSEK_ROOT_DIR}/h/")
find_library (MOSEK_LIBRARY1  libfusion64.a   HINTS "${MOSEK_ROOT_DIR}/bin/")
find_library (MOSEK_LIBRARY2  libmosek64.a   HINTS "${MOSEK_ROOT_DIR}/bin/")


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MOSEK DEFAULT_MSG MOSEK_LIBRARY1 MOSEK_LIBRARY2 MOSEK_INCLUDE_DIR)

if(MOSEK_FOUND)
	message("—- Found Sdpa under ${MOSEK_INCLUDE_DIR}")
    set(MOSEK_INCLUDE_DIRS ${MOSEK_INCLUDE_DIR})
    set(MOSEK_LIBRARIES  ${MOSEK_LIBRARY1} ${MOSEK_LIBRARY2})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(MOSEK_LIBRARIES "${MOSEK_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(MOSEK_FOUND)

mark_as_advanced(MOSEK_LIBRARY MOSEK_INCLUDE_DIR)
