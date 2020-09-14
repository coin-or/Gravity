# Search supplied hint directories first if supplied.
set(GoICP_ROOT_DIR "$ENV{GoICP_ROOT_DIR}" CACHE PATH "GoICP root directory.")
if(GoICP_ROOT_DIR)
 message("Looking for GoICP in ${GoICP_ROOT_DIR}")
else(GoICP_ROOT_DIR)
 message("GoICP_ROOT_DIR not provided.")
endif(GoICP_ROOT_DIR)

find_path(GoICP_INCLUDE_DIR
  NAMES jly_goicp.h
  HINTS /usr/local/include
  HINTS /usr/include
  HINTS ${THIRDPARTY_INSTALL_PATH}/GoICP/include
  HINTS ${THIRDPARTY_INSTALL_PATH}/include)

# Set standard CMake FindPackage variables if found.
if (GoICP_FOUND)
  set(GoICP_INCLUDE_DIRS ${GoICP_INCLUDE_DIR})
else (GoICP_FOUND)
 message("Cannot find GoICP, will try pulling it from github.")
endif (GoICP_FOUND)
