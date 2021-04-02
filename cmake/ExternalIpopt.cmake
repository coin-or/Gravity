# Create download URL derived from version number.
set(IPOPT_HOME https://www.coin-or.org/Tarballs/Ipopt)
set(IPOPT_DOWNLOAD_URL ${IPOPT_HOME}/Ipopt-3.12.13.tgz)
unset(IPOPT_HOME)

# Download and build the IPOPT library and add its properties to the third party arguments.
set(IPOPT_ROOT_DIR ${THIRDPARTY_INSTALL_PATH}/Install/ipopt/build CACHE INTERNAL "")
if(WIN32)
ExternalProject_Add(ipopt
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND curl -k -L https://github.com/IDAES/idaes-ext/releases/download/2.4.1/idaes-solvers-windows-64.tar.gz -o idaes-solvers-windows-64.tar.gz && tar -xvzf idaes-solvers-windows-64.tar.gz -C ${PROJECT_SOURCE_DIR}/bin/Release
    URL ${IPOPT_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${IPOPT_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
else()
ExternalProject_Add(ipopt
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && curl -k -L ${IPOPT_DOWNLOAD_URL} -o Ipopt.tar.gz && tar -xzf Ipopt.tar.gz && rm -fr ./Install/ipopt && mv Ipopt-3.12.13 ./Install/ipopt && cd ./Install/ipopt && mkdir build && cd ./ThirdParty/Mumps && export HTTP_PROXY=$ENV{HTTP_PROXY} && ./get.Mumps && cd ../../build && ../configure --prefix=${IPOPT_ROOT_DIR} && make install -j24
    URL ${IPOPT_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${IPOPT_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
endif()
list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DIPOPT_ROOT_DIR:PATH=${IPOPT_ROOT_DIR}")
set(IPOPT_INCLUDE_DIRS ${THIRDPARTY_INSTALL_PATH}/Install/ipopt/build/include/coin)
include_directories(${IPOPT_INCLUDE_DIRS})
if(APPLE)
find_library(IPOPT_LIBRARIES
        libipopt.dylib
        HINTS /usr/local/lib
        HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
        HINTS ${IPOPT_ROOT_DIR}/lib
)
elseif(UNIX)
find_library(IPOPT_LIBRARIES
        libipopt.so
        HINTS /usr/local/lib
        HINTS ${IPOPT_ROOT_DIR}/lib
        HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
)
endif()
set(LIBS ${LIBS} ${IPOPT_LIBRARIES})
unset(IPOPT_DOWNLOAD_URL)
unset(IPOPT_ROOT)
