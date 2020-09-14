# Create download URL derived from version number.
set(GoICP_HOME https://github.com/yangjiaolong/Go-ICP/archive)
set(GoICP_DOWNLOAD_URL ${GoICP_HOME}/master.tar.gz)
unset(GoICP_HOME)

# Download and build the GoICP library and add its properties to the third party arguments.
set(GoICP_ROOT ${THIRDPARTY_INSTALL_PATH} CACHE INTERNAL "")
if(APPLE)
ExternalProject_Add(GoICP
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && curl -k -L ${GoICP_DOWNLOAD_URL} -o GoICP.tar.gz && tar -xzf GoICP.tar.gz && mv Go-ICP-master GoICP && rm -fr ./Install/GoICP && sed -i.bu 13d GoICP/jly_3ddt.cpp && mv GoICP ./Install && cd Install/GoICP && cmake . && make GoICP
    URL ${GoICP_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${GoICP_ROOT}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
elseif(UNIX)
ExternalProject_Add(GoICP
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && curl -k -L ${GoICP_DOWNLOAD_URL} -o GoICP.tar.gz && tar -xzf GoICP.tar.gz && mv Go-ICP-master GoICP && rm -fr ./Install/GoICP && sed -i '13d' jly_3ddt.cpp && mv GoICP ./Install && cd Install/GoICP && cmake . && make GoICP
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${GoICP_ROOT}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
endif()

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DGoICP_ROOT:PATH=${GoICP_ROOT}")
unset(GoICP_DOWNLOAD_URL)
unset(GoICP_ROOT)
