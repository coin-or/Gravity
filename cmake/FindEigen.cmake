# Search supplied hint directories first if supplied.
set(EIGEN_ROOT_DIR "$ENV{EIGEN_ROOT_DIR}" CACHE PATH "Eigen root directory.")
if(EIGEN_ROOT_DIR)
 message("Looking for Eigen in ${EIGEN_ROOT_DIR}")
else(EIGEN_ROOT_DIR)
 message("EIGEN_ROOT_DIR not provided.")
endif(EIGEN_ROOT_DIR)

find_path(EIGEN_INCLUDE_DIR
  NAMES Eigen/Core
  HINTS /usr/local/include
  HINTS /usr/local/include/eigen3
  HINTS /usr/include
  HINTS /usr/include/eigen3
  HINTS ${THIRDPARTY_INSTALL_PATH}/include)

# Extract Eigen version from Eigen/src/Core/util/Macros.h
if (EIGEN_INCLUDE_DIR)
  set(EIGEN_VERSION_FILE ${EIGEN_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h)
  if (NOT EXISTS ${EIGEN_VERSION_FILE})
    eigen_report_not_found(
      "Could not find file: ${EIGEN_VERSION_FILE} "
      "containing version information in Eigen install located at: "
      "${EIGEN_INCLUDE_DIR}.")
  else (NOT EXISTS ${EIGEN_VERSION_FILE})
    file(READ ${EIGEN_VERSION_FILE} EIGEN_VERSION_FILE_CONTENTS)

    string(REGEX MATCH "#define EIGEN_WORLD_VERSION [0-9]+"
      EIGEN_WORLD_VERSION "${EIGEN_VERSION_FILE_CONTENTS}")
    string(REGEX REPLACE "#define EIGEN_WORLD_VERSION ([0-9]+)" "\\1"
      EIGEN_WORLD_VERSION "${EIGEN_WORLD_VERSION}")

    string(REGEX MATCH "#define EIGEN_MAJOR_VERSION [0-9]+"
      EIGEN_MAJOR_VERSION "${EIGEN_VERSION_FILE_CONTENTS}")
    string(REGEX REPLACE "#define EIGEN_MAJOR_VERSION ([0-9]+)" "\\1"
      EIGEN_MAJOR_VERSION "${EIGEN_MAJOR_VERSION}")

    string(REGEX MATCH "#define EIGEN_MINOR_VERSION [0-9]+"
      EIGEN_MINOR_VERSION "${EIGEN_VERSION_FILE_CONTENTS}")
    string(REGEX REPLACE "#define EIGEN_MINOR_VERSION ([0-9]+)" "\\1"
      EIGEN_MINOR_VERSION "${EIGEN_MINOR_VERSION}")

    # This is on a single line s/t CMake does not interpret it as a list of
    # elements and insert ';' separators which would result in 3.;2.;0 nonsense.
    set(EIGEN_VERSION "${EIGEN_WORLD_VERSION}.${EIGEN_MAJOR_VERSION}.${EIGEN_MINOR_VERSION}")
  endif (NOT EXISTS ${EIGEN_VERSION_FILE})
endif (EIGEN_INCLUDE_DIR)

# Set standard CMake FindPackage variables if found.
if (EIGEN_FOUND)
  set(EIGEN_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR})
else (EIGEN_FOUND)
 message("Cannot find Eigen, will try pulling it from github.")
endif (EIGEN_FOUND)
