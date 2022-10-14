if(CMAKE_XCODE_BUILD_SYSTEM VERSION_GREATER_EQUAL 12)
  cmake_policy(SET CMP0114 NEW)
else()
  cmake_policy(SET CMP0114 OLD)
endif()

set(HiGHS_DOWNLOAD_URL https://github.com/ERGO-Code/HiGHS)



# Download and build the HiGHS library and add its properties to the third party arguments.
set(HiGHS_ROOT_DIR ${THIRDPARTY_INSTALL_PATH}/Install/HiGHS CACHE INTERNAL "")
ExternalProject_Add(HiGHS
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && git clone ${HiGHS_DOWNLOAD_URL} && rm -fr ./Install/HiGHS && mv HiGHS ./Install/HiGHS && cd ./Install/HiGHS && mkdir build && cd build && cmake .. && make -j
    URL ${HiGHS_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${HiGHS_ROOT_DIR}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DHiGHS_ROOT_DIR:PATH=${HiGHS_ROOT_DIR}")
set(HiGHS_INCLUDE_DIRS ${THIRDPARTY_INSTALL_PATH}/Install/HiGHS/src ${THIRDPARTY_INSTALL_PATH}/Install/HiGHS/build)
include_directories(${HiGHS_INCLUDE_DIRS})
set(LIBS ${LIBS} ${HiGHS_LIBRARIES})
