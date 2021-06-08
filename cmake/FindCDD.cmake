set(CDD_ROOT_DIR "$ENV{CDD_ROOT_DIR}" CACHE PATH "CDD root directory.")
message("Looking for CDD in ${CDD_ROOT_DIR}")

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DCDD_ROOT_DIR:PATH=${CDD_ROOT_DIR}")
include_directories(${CDD_INCLUDE_DIRS})

set(CDD_INCLUDE_DIRS ${CDD_ROOT_DIR}/lib-src)

if(APPLE)
find_library(CDD_LIBRARY 
	libcdd.0.dylib
	HINTS /usr/local/lib
	HINTS ${CDD_ROOT_DIR}/lib/lib
	message("Looking for CDDLib in ${CDD_ROOT_DIR}/lib/lib")
)
elseif(UNIX)
find_library(CDD_LIBRARY 
	libcdd.so.0
	HINTS /usr/local/lib
	HINTS ${CDD_ROOT_DIR}/lib/lib
)
endif()

message("CDD_Library is ${CDD_LIBRARY}")

set(CDD_LIBRARIES ${CDD_LIBRARY})
set(LIBS ${LIBS} ${CDD_LIBRARIES})

