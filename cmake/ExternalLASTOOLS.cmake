# Create download URL derived from version number.
set(LASTOOLS_DOWNLOAD_URL https://github.com/LAStools/LAStools/archive/refs/heads/master.zip)

# Download and build the Lastools library and add its properties to the third party arguments.
set(LASTOOLS_ROOT_DIR ${THIRDPARTY_INSTALL_PATH}/Install/lastools CACHE INTERNAL "")
ExternalProject_Add(lastools
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTP_PROXY=$ENV{HTTP_PROXY} && export HTTPS_PROXY=$ENV{HTTPS_PROXY} && export http_proxy=$ENV{HTTP_PROXY} && curl -k -L ${LASTOOLS_DOWNLOAD_URL} -o Lastools.zip && unzip Lastools.zip && rm -fr ./Install/lastools && mv LAStools-master ./Install/lastools && cd ./Install/lastools && make -j24
    URL ${LASTOOLS_DOWNLOAD_URL}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DLASTOOLS_ROOT_DIR:PATH=${LASTOOLS_ROOT_DIR}")
set(LASlib_INCLUDE_DIRS ${LASTOOLS_ROOT_DIR}/LASlib/inc)
include_directories(${LASlib_INCLUDE_DIRS})
find_library(LASlib_LIBRARY1
        liblas.a
        HINTS /usr/local/lib
        HINTS ${LASTOOLS_ROOT_DIR}/LASlib/lib
)


set(LASlib_LIBRARIES  ${LASlib_LIBRARY1})
set(LIBS ${LIBS} ${LASlib_LIBRARIES})
unset(LASTOOLS_DOWNLOAD_URL)
unset(LASTOOLS_ROOT)
