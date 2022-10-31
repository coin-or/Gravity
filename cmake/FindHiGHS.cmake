set(HiGHS_ROOT_DIR "$ENV{HiGHS_ROOT_DIR}" CACHE PATH "HiGHS root directory.")
message("Looking for HiGHS in ${HiGHS_ROOT_DIR}")


find_path(HiGHS_INCLUDE_DIR
	NAMES Highs.h	
	HINTS ${HiGHS_ROOT_DIR}/include/HiGHS
	HINTS ${PROJECT_SOURCE_DIR}/third_party/HiGHS/build/include
	HINTS /usr/local/include/
)

if(WIN32)

	find_library(HiGHS_LIBRARY 
		libhighs.lib
		HINTS /usr/local/lib
		HINTS ${PROJECT_SOURCE_DIR}/third_party/HiGHS/build/lib
		HINTS ${HiGHS_ROOT_DIR}/build/lib
	)
elseif(APPLE)
	find_library(HiGHS_LIBRARY 
		libhighs.dylib
		HINTS /usr/local/lib
		HINTS ${PROJECT_SOURCE_DIR}/third_party/HiGHS/build/lib
		HINTS ${HiGHS_ROOT_DIR}/build/lib
	)
elseif(UNIX)
	find_library(HiGHS_LIBRARY 
		libhighs.so
		HINTS /usr/local/lib
		HINTS ${PROJECT_SOURCE_DIR}/third_party/HiGHS/build/lib
		HINTS ${HiGHS_ROOT_DIR}/build/lib
	)
endif(WIN32)
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HiGHS DEFAULT_MSG HiGHS_LIBRARY HiGHS_INCLUDE_DIR)

if(HiGHS_FOUND)
	message("—- Found HiGHS under ${HiGHS_INCLUDE_DIR}")
    set(HiGHS_INCLUDE_DIRS ${HiGHS_INCLUDE_DIR})
    set(HiGHS_LIBRARIES ${HiGHS_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(HiGHS_LIBRARIES "${HiGHS_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
else (HiGHS_FOUND)
 message("Cannot find HiGHS, will try pulling it from github.")
endif(HiGHS_FOUND)

mark_as_advanced(HiGHS_LIBRARY HiGHS_INCLUDE_DIR)
