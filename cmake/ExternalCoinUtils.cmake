# Create download URL derived from version number.
set(CoinUtils_DOWNLOAD_URL https://github.com/coin-or/CoinUtils)
unset(CoinUtils_HOME)

# Download and build the CoinUtils library and add its properties to the third party arguments.
set(CoinUtils_ROOT_DIR ${THIRDPARTY_INSTALL_PATH}/Install/CoinUtils/build CACHE INTERNAL "")
if(WIN32)
ExternalProject_Add(coinutils
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND git clone ${CoinUtils_DOWNLOAD_URL} && rmdir -fr ./Install/CoinUtils && move CoinUtils ./Install/CoinUtils && cd ./Install/CoinUtils && mkdir build && cd build && ../configure --prefix=${CoinUtils_ROOT_DIR} && make -j24 && make install
    URL ${CoinUtils_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CoinUtils_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
else()
ExternalProject_Add(coinutils
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && git clone ${CoinUtils_DOWNLOAD_URL} && rm -fr ./Install/CoinUtils && mv CoinUtils ./Install/CoinUtils && cd ./Install/CoinUtils && mkdir build && cd build && ../configure --prefix=${CoinUtils_ROOT_DIR} && make -j24 && make install
    URL ${CoinUtils_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CoinUtils_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
endif()
list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DCoinUtils_ROOT_DIR:PATH=${CoinUtils_ROOT_DIR}")
set(CoinUtils_INCLUDE_DIRS ${THIRDPARTY_INSTALL_PATH}/Install/CoinUtils/build/include/coin-or)
include_directories(${CoinUtils_INCLUDE_DIRS})
if(APPLE)
find_library(CoinUtils_LIBRARY
        libCoinUtils.dylib
        HINTS ${CoinUtils_ROOT_DIR}/lib
)
elseif(UNIX)
find_library(CoinUtils_LIBRARY
        libCoinUtils.so
        HINTS ${CoinUtils_ROOT_DIR}/lib
)
endif()
set(LIBS ${LIBS} ${CoinUtils_LIBRARY})
unset(CoinUtils_DOWNLOAD_URL)
unset(CoinUtils_ROOT)
