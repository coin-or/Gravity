# Create download URL derived from version number.
set(IPOPT_HOME https://www.coin-or.org/Tarballs/Ipopt)
set(IPOPT_DOWNLOAD_URL ${IPOPT_HOME}/Ipopt-3.12.13.tgz)
unset(IPOPT_HOME)

set(IPOPT_ROOT_DIR ${PROJECT_SOURCE_DIR}/thirdparty/Ipopt CACHE INTERNAL "")

if(WIN32)
ExternalProject_Add(ipopt
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND curl -k -L https://github.com/IDAES/idaes-ext/releases/download/2.4.1/idaes-solvers-windows-64.tar.gz -o idaes-solvers-windows-64.tar.gz && tar -xvzf idaes-solvers-windows-64.tar.gz -C ${IPOPT_ROOT_DIR}
    URL ${IPOPT_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${IPOPT_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
add_custom_command(
  TARGET ipopt POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy 
  ${IPOPT_ROOT_DIR}/libipopt-3.dll ${PROJECT_SOURCE_DIR}/bin/Release/libipopt-3.dll)

add_custom_command(
TARGET ipopt POST_BUILD
COMMAND ${CMAKE_COMMAND} -E copy 
${IPOPT_ROOT_DIR}/libblas.dll ${PROJECT_SOURCE_DIR}/bin/Release/libblas.dll)

add_custom_command(
TARGET ipopt POST_BUILD
COMMAND ${CMAKE_COMMAND} -E copy 
${IPOPT_ROOT_DIR}/liblapack.dll ${PROJECT_SOURCE_DIR}/bin/Release/liblapack.dll)

add_custom_command(
TARGET ipopt POST_BUILD
COMMAND ${CMAKE_COMMAND} -E copy 
${IPOPT_ROOT_DIR}/libgfortran-5.dll ${PROJECT_SOURCE_DIR}/bin/Release/libgfortran-5.dll)

else()

# Download and build the IPOPT library and add its properties to the third party arguments.
ExternalProject_Add(ipopt
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && curl -k -L ${IPOPT_DOWNLOAD_URL} -o Ipopt.tar.gz && tar -xzf Ipopt.tar.gz && mv Ipopt-3.12.13 ${IPOPT_ROOT_DIR} && cd ${IPOPT_ROOT_DIR} && mkdir build && cd ./ThirdParty/Mumps && export HTTP_PROXY=$ENV{HTTP_PROXY} && ./get.Mumps && cd ../../build && ../configure --prefix=${IPOPT_ROOT_DIR} && make install -j24
    URL ${IPOPT_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${IPOPT_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
endif()
list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DIPOPT_ROOT_DIR:PATH=${IPOPT_ROOT_DIR}")
set(IPOPT_INCLUDE_DIRS ${IPOPT_ROOT_DIR}/include/coin)
include_directories(${IPOPT_INCLUDE_DIRS})
if(APPLE)
find_library(IPOPT_LIBRARIES
        libipopt.dylib
        HINTS /usr/local/lib
        HINTS ${IPOPT_ROOT_DIR}/build/lib
        HINTS ${IPOPT_ROOT_DIR}/lib
        HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
)
elseif(UNIX)
find_library(IPOPT_LIBRARIES
        libipopt.so
        HINTS /usr/local/lib
        HINTS ${IPOPT_ROOT_DIR}/build/lib
        HINTS ${IPOPT_ROOT_DIR}/lib
        HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinIpopt/build/lib
)
endif()
set(LIBS ${LIBS} ${IPOPT_LIBRARIES})
unset(IPOPT_DOWNLOAD_URL)
unset(IPOPT_ROOT)
