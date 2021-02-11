# Create download URL derived from version number.
set(QHULL_DOWNLOAD_URL https://github.com/qhull/qhull/archive/master.zip)

# Download and build the Qhull library and add its properties to the third party arguments.
set(QHULL_ROOT_DIR ${THIRDPARTY_INSTALL_PATH}/Install/qhull/build CACHE INTERNAL "")
ExternalProject_Add(qhull
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && curl -k -L ${QHULL_DOWNLOAD_URL} -o Qhull.zip && unzip Qhull.zip && rm -fr ./Install/qhull && mv qhull-master ./Install/qhull && cd ./Install/qhull/build && cmake -B ${QHULL_ROOT_DIR} .. && make -j24
    URL ${QHULL_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${QHULL_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DQHULL_ROOT_DIR:PATH=${QHULL_ROOT_DIR}")
set(QHULL_INCLUDE_DIRS ${THIRDPARTY_INSTALL_PATH}/Install/qhull/src)
include_directories(${QHULL_INCLUDE_DIRS})
find_library(QHULL_LIBRARY1
        libqhullcpp.a
        HINTS /usr/local/lib
        HINTS ${QHULL_ROOT_DIR}
)

find_library(QHULL_LIBRARY2
        libqhullstatic.a
        HINTS /usr/local/lib
        HINTS ${QHULL_ROOT_DIR}
)

find_library(QHULL_LIBRARY3
        libqhullstatic_r.a
        HINTS /usr/local/lib
        HINTS ${QHULL_ROOT_DIR}
)
#set(QHULL_LIBRARIES  ${QHULL_LIBRARY1} ${QHULL_LIBRARY2} ${QHULL_LIBRARY3})
set(QHULL_LIBRARIES  ${QHULL_LIBRARY1} ${QHULL_LIBRARY3})
set(LIBS ${LIBS} ${QHULL_LIBRARIES})
unset(QHULL_DOWNLOAD_URL)
unset(QHULL_ROOT)
