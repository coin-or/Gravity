# Create download URL derived from version number.
set(LASTOOLS_DOWNLOAD_URL https://github.com/LAStools/LAStools.git)


# Download and build the LASTOOLS library and add its properties to the third party arguments.
set(LASTOOLS_ROOT_DIR ${THIRDPARTY_INSTALL_PATH}/Install/LASTOOLS CACHE INTERNAL "")

if(WIN32)
# On Windows, download the Ipopt dlls from https://github.com/IDAES and the LASTOOLS dlls from https://github.com/coin-or/LASTOOLS
ExternalProject_Add(LASTOOLS
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND curl -k -L LAStools.tar.gz -o LAStools.tar.gz && tar -xvzf LAStools.tar.gz -C ${LASTOOLS_ROOT_DIR}
    URL ${LASTOOLS_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${LASTOOLS_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
else()
ExternalProject_Add(LASTOOLS
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && git clone ${LASTOOLS_DOWNLOAD_URL} && rm -fr ./Install/LASTOOLS && mv LASTOOLS ./Install/LASTOOLS && cd ./Install/LASTOOLS && make
    URL ${MP_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${MP_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
endif()

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DLASTOOLS_ROOT_DIR:PATH=${LASTOOLS_ROOT_DIR}")
set(LASLIB_INCLUDE_DIRS ${THIRDPARTY_INSTALL_PATH}/Install/LASTOOLS/LASlib/inc ${THIRDPARTY_INSTALL_PATH}/Install/LASTOOLS/LASzip/src)
set(LASLIB_LIBRARY ${THIRDPARTY_INSTALL_PATH}/Install/LASTOOLS/LASlib/lib)
include_directories(${LASLIB_INCLUDE_DIRS})
set(LIBS ${LIBS} ${LASLIB_LIBRARY})
