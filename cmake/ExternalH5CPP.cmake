if(CMAKE_XCODE_BUILD_SYSTEM VERSION_GREATER_EQUAL 12)
  cmake_policy(SET CMP0114 NEW)
else()
  cmake_policy(SET CMP0114 OLD)
endif()

set(H5CPP_DOWNLOAD_URL https://github.com/ess-dmsc/h5cpp.git)



# Download and build the H5CPP library and add its properties to the third party arguments.
set(H5CPP_ROOT_DIR ${THIRDPARTY_INSTALL_PATH}/Install/H5CPP/build CACHE INTERNAL "")
ExternalProject_Add(h5cpp
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && git clone ${H5CPP_DOWNLOAD_URL} && rm -fr ./Install/H5CPP && mv H5CPP ./Install/H5CPP && cd ./Install/H5CPP && mkdir build && cd build && cmake .. -DH5CPP_WITH_BOOST=OFF -DH5CPP_CONAN=DISABLE -DH5CPP_DISABLE_TESTS=ON && make
    URL ${H5CPP_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${H5CPP_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DH5CPP_ROOT_DIR:PATH=${H5CPP_ROOT_DIR}")
set(H5CPP_INCLUDE_DIRS ${THIRDPARTY_INSTALL_PATH}/Install/H5CPP/build/src)
include_directories(${H5CPP_INCLUDE_DIRS})
if(APPLE)
find_library(H5CPP_LIBRARIES
        libh5cpp.dylib
        HINTS /usr/local/lib
        HINTS ${PROJECT_SOURCE_DIR}/third_party/Install/H5CPP/build/lib
        HINTS ${H5CPP_ROOT_DIR}/lib
)
#elseif(UNIX)
#find_library(H5CPP_LIBRARIES
#        libH5CPP.so
#        HINTS /usr/local/lib
#        HINTS ${H5CPP_ROOT_DIR}/lib
#        HINTS ${PROJECT_SOURCE_DIR}/third_party/CoinH5CPP/build/lib
#)
endif()
set(LIBS ${LIBS} ${H5CPP_LIBRARIES})
