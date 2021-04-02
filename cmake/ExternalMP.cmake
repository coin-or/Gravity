# Create download URL derived from version number.
set(MP_DOWNLOAD_URL https://github.com/ampl/mp.git)
unset(MP_HOME)

# Download and build the MP library and add its properties to the third party arguments.
set(MP_ROOT_DIR ${THIRDPARTY_INSTALL_PATH}/Install/MP CACHE INTERNAL "")
if(WIN32)
ExternalProject_Add(mp
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND git clone ${MP_DOWNLOAD_URL} && rmdir -fr ./Install/MP && move mp ./Install/MP && cd ./Install/MP && git submodule init && git submodule update && mkdir build && cd build && cmake -DBUILD=all .. && make -j24
    URL ${MP_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${MP_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
else()
ExternalProject_Add(mp
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && git clone ${MP_DOWNLOAD_URL} && rm -fr ./Install/MP && mv mp ./Install/MP && cd ./Install/MP && git submodule init && git submodule update && mkdir build && cd build && cmake -DBUILD=all .. && make -j24
    URL ${MP_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${MP_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
endif()

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DMP_ROOT_DIR:PATH=${MP_ROOT_DIR}")
set(MP_INCLUDE_DIRS ${THIRDPARTY_INSTALL_PATH}/Install/MP/include)
include_directories(${MP_INCLUDE_DIRS})
find_library(MP_LIBRARY
        libmp.a
        HINTS ${MP_ROOT_DIR}/build/lib
)
set(LIBS ${LIBS} ${MP_LIBRARY})
unset(MP_DOWNLOAD_URL)
unset(MP_ROOT)
