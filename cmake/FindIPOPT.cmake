set(IPOPT_ROOT_DIR "$ENV{IPOPT_ROOT_DIR}" CACHE PATH "IPOPT root directory.")
message("Looking for Ipopt in ${IPOPT_ROOT_DIR}")


include(FindPackageHandleStandardArgs)

if(WIN32)

find_path(IPOPT_INCLUDE_DIR
	NAMES IpNLP.hpp 
	HINTS /usr/local/include/coin
	HINTS ${IPOPT_ROOT_DIR}/include/coin
	HINTS ${IPOPT_ROOT_DIR}/include
	HINTS ${PROJECT_SOURCE_DIR}/thirdparty/Ipopt/include/coin-or
)

find_library(IPOPT_LIBRARY 
	libipopt-3.lib
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/thirdparty/Ipopt
	HINTS ${IPOPT_ROOT_DIR}/lib
)

find_library(IPOPT_LIBRARY2
	# libipopt-3.dll
	libpynumero_ASL.dll.a
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/thirdparty/Ipopt
	HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
	HINTS ${IPOPT_ROOT_DIR}/lib
) 

find_package_handle_standard_args(IPOPT DEFAULT_MSG IPOPT_LIBRARY IPOPT_LIBRARY2 IPOPT_INCLUDE_DIR)
elseif(APPLE)
find_path(IPOPT_INCLUDE_DIR
	NAMES IpNLP.hpp 
	HINTS /usr/local/include/coin
	HINTS ${IPOPT_ROOT_DIR}/include/coin
	HINTS ${IPOPT_ROOT_DIR}/include
)
find_library(IPOPT_LIBRARY 
	libipopt.dylib
	HINTS /usr/local/lib
	HINTS ${IPOPT_ROOT_DIR}/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
)
find_package_handle_standard_args(IPOPT DEFAULT_MSG IPOPT_LIBRARY IPOPT_INCLUDE_DIR)
elseif(UNIX)
find_path(IPOPT_INCLUDE_DIR
	NAMES IpNLP.hpp 
	HINTS /usr/local/include/coin
	HINTS ${IPOPT_ROOT_DIR}/include/coin
	HINTS ${IPOPT_ROOT_DIR}/include
	HINTS ${PROJECT_SOURCE_DIR}/thirdparty/Ipopt/include/coin-or
)
find_library(IPOPT_LIBRARY 
	libipopt.so
	HINTS /usr/local/lib
	HINTS ${IPOPT_ROOT_DIR}/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
	HINTS ${PROJECT_SOURCE_DIR}/thirdparty/Ipopt
)
find_package_handle_standard_args(IPOPT DEFAULT_MSG IPOPT_LIBRARY IPOPT_INCLUDE_DIR)
endif(WIN32)

if(IPOPT_FOUND)
	message("—- Found Ipopt include dir under ${IPOPT_INCLUDE_DIR}")
	message("—- Found Ipopt lib at ${IPOPT_LIBRARY}")
    set(IPOPT_INCLUDE_DIRS ${IPOPT_INCLUDE_DIR})
    set(IPOPT_LIBRARIES ${IPOPT_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(IPOPT_LIBRARIES "${IPOPT_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
else (IPOPT_FOUND)
 message("Cannot find Ipopt, will try pulling it from github.")
endif(IPOPT_FOUND)

mark_as_advanced(IPOPT_LIBRARY IPOPT_INCLUDE_DIR)

