if(CMAKE_XCODE_BUILD_SYSTEM VERSION_GREATER_EQUAL 12)
  cmake_policy(SET CMP0114 NEW)
else()
  cmake_policy(SET CMP0114 OLD)
endif()

# Create download URL derived from version number.
set(SMTLIB_DOWNLOAD_URL https://es-static.fbk.eu/people/griggio/misc/downloads/smtlib2parser-1.4.tar.gz)
unset(SMTLIB_HOME)

# Download and build the SMTLIB library and add its properties to the third party arguments.
set(SMTLIB_ROOT_DIR ${THIRDPARTY_INSTALL_PATH}/Install/SMTLIB CACHE INTERNAL "")
if(WIN32)
ExternalProject_Add(SMTLIB
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND curl -k -L ${SMTLIB_DOWNLOAD_URL} -o SMTLIB.tar.gz && tar -xzf SMTLIB.tar.gz && move smtlib2parser-1.4 SMTLIB && rmdir -fr ./Install/SMTLIB && move SMTLIB ./Install/SMTLIB && cd ./Install/SMTLIB && make libsmtlib2parser.a
    URL ${MP_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${MP_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
else()
ExternalProject_Add(SMTLIB
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && curl -k -L ${SMTLIB_DOWNLOAD_URL} -o SMTLIB.tar.gz && tar -xzf SMTLIB.tar.gz && mv smtlib2parser-1.4 SMTLIB && rm -fr ./Install/SMTLIB && mv SMTLIB ./Install/SMTLIB && cd ./Install/SMTLIB && mv Makefile.in Makefile && make libsmtlib2parser.a
    URL ${SMTLIB_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${SMTLIB_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
endif()

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DSMTLIB_ROOT_DIR:PATH=${SMTLIB_ROOT_DIR}")
set(SMTLIB_INCLUDE_DIRS ${THIRDPARTY_INSTALL_PATH}/Install/SMTLIB)
include_directories(${SMTLIB_INCLUDE_DIRS})
find_library(SMTLIB_LIBRARY
        libsmtlib2parser.a
        HINTS ${SMTLIB_ROOT_DIR}
)
set(LIBS ${LIBS} ${SMTLIB_LIBRARY})
unset(SMTLIB_DOWNLOAD_URL)
unset(SMTLIB_ROOT)
