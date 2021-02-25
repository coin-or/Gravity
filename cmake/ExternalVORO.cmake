# Create download URL derived from version number.
set(VORO_DOWNLOAD_URL http://math.lbl.gov/voro++/download/dir/voro++-0.4.6.tar.gz)

# Download and build the Voro library and add its properties to the third party arguments.
set(VORO_ROOT_DIR ${THIRDPARTY_INSTALL_PATH}/Install/voro CACHE INTERNAL "")
ExternalProject_Add(voro
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && curl -k -L ${VORO_DOWNLOAD_URL} -o Voro.tar.gz && tar -xvf Voro.tar.gz && rm -fr ./Install/voro && mv voro++-0.4.6 ./Install/voro && cd ./Install/voro/build && make -j24
    URL ${VORO_DOWNLOAD_URL}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DVORO_ROOT_DIR:PATH=${VORO_ROOT_DIR}")
set(VORO_INCLUDE_DIRS ${VORO_ROOT_DIR}/src)
include_directories(${VORO_INCLUDE_DIRS})
find_library(VORO_LIBRARY1
        libvoro++.a
        HINTS /usr/local/lib
        HINTS ${VORO_ROOT_DIR}/src
)


set(VORO_LIBRARIES  ${VORO_LIBRARY1})
set(LIBS ${LIBS} ${VORO_LIBRARIES})
unset(VORO_DOWNLOAD_URL)
unset(VORO_ROOT)
