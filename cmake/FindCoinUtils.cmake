set(CoinUtils_ROOT_DIR "$ENV{CoinUtils_ROOT_DIR}" CACHE PATH "CoinUtils root directory.")
message("Looking for CoinUtils in ${CoinUtils_ROOT_DIR}")


find_path(CoinUtils_INCLUDE_DIR
	NAMES CoinMpsIO.hpp
	HINTS /usr/local/include
	HINTS ${CoinUtils_ROOT_DIR}/include/coin-or
	HINTS ${CoinUtils_ROOT_DIR}/include/coin
	HINTS ${CoinUtils_ROOT_DIR}/include/coinr
	HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinUtils/build/include/coin-or
)

if(APPLE)
find_library(CoinUtils_LIBRARY 
	libCoinUtils.dylib
	HINTS /usr/local/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinUtils/build/lib
	HINTS ${CoinUtils_ROOT_DIR}/lib
)
elseif(UNIX)
find_library(CoinUtils_LIBRARY 
	libCoinUtils.so
	HINTS /usr/local/lib
	HINTS ${CoinUtils_ROOT_DIR}/lib
	HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinUtils/build/lib
)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CoinUtils DEFAULT_MSG CoinUtils_LIBRARY CoinUtils_INCLUDE_DIR)

if(CoinUtils_FOUND)
	message("—- Found CoinUtils under ${CoinUtils_INCLUDE_DIR}")
    set(CoinUtils_INCLUDE_DIRS ${CoinUtils_INCLUDE_DIR})
    set(CoinUtils_LIBRARIES ${CoinUtils_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(CoinUtils_LIBRARIES "${CoinUtils_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
else (CoinUtils_FOUND)
 message("Cannot find CoinUtils, will try pulling it from github.")
endif(CoinUtils_FOUND)

mark_as_advanced(CoinUtils_LIBRARY CoinUtils_INCLUDE_DIR)
