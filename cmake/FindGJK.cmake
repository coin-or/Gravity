set(GJK_ROOT_DIR "$ENV{GJK_ROOT_DIR}" CACHE PATH "GJK root directory.")
message("Looking for GJK in ${GJK_ROOT_DIR}")

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DGJK_ROOT_DIR:PATH=${GJK_ROOT_DIR}")
include_directories(${GJK_INCLUDE_DIRS})

find_path(GJK_INCLUDE_DIRS
	NAMES openGJK.h
	HINTS ${GJK_ROOT_DIR}/lib/include/openGJK
	HINTS ${GJK_ROOT_DIR}/include/openGJK
	HINTS /usr/local/include/
)

find_library(GJK_LIBRARY 
	libopenGJKlib.a
	HINTS /usr/local/lib
	HINTS ${GJK_ROOT_DIR}/build
	HINTS ${GJK_ROOT_DIR}/build/lib
	HINTS ${GJK_ROOT_DIR}/lib
)

if(APPLE)
if(NOT GJK_LIBRARY)
find_library(GJK_LIBRARY 
	libopenGJKlib.dylib
	HINTS /usr/local/lib
	HINTS ${GJK_ROOT_DIR}/build
	HINTS ${GJK_ROOT_DIR}/build/lib
	HINTS ${GJK_ROOT_DIR}/lib
)
elseif(UNIX)
find_library(GJK_LIBRARY 
	libopenGJKlib.so
	HINTS /usr/local/lib
	HINTS ${GJK_ROOT_DIR}/build
	HINTS ${GJK_ROOT_DIR}/build/lib
	HINTS ${GJK_ROOT_DIR}/lib
)
endif()
endif()
message("GJK_Library is ${GJK_LIBRARY}")

set(GJK_LIBRARIES ${GJK_LIBRARY})
set(LIBS ${LIBS} ${GJK_LIBRARIES})

