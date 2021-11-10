set(CLP_ROOT_DIR "$ENV{CLP_ROOT_DIR}" CACHE PATH "CLP root directory.")
message("Looking for CLP in ${CLP_ROOT_DIR} folder")

find_path(CLP_INCLUDE_DIR  
    NAMES ClpSimplex.hpp 
    HINTS ${CLP_ROOT_DIR}/include/coin/
    )

find_library(CLP_LIBRARY
    NAMES libClp.dylib
    HINTS ${CLP_ROOT_DIR}/lib/
    )
find_library(UTILS_LIBRARY
  NAMES libCoinUtils.dylib
  HINTS ${CLP_ROOT_DIR}/lib/
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CLP DEFAULT_MSG UTILS_LIBRARY CLP_LIBRARY CLP_INCLUDE_DIR)

if(CLP_FOUND)
    message("—- Found CLP under ${CLP_INCLUDE_DIR}")
    set(CLP_INCLUDE_DIRS ${CLP_INCLUDE_DIR})
    set(CLP_LIBRARIES ${CLP_LIBRARY} ${UTILS_LIBRARY} ${CLP_ROOT_DIR}/lib)
    link_directories(${CLP_ROOT_DIR}/lib/)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(CLP_LIBRARIES "${CLP_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(CLP_FOUND)

mark_as_advanced(CLP_LIBRARY UTILS_LIBRARY CLP_INCLUDE_DIR)
