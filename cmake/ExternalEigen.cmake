# Create download URL derived from version number.
set(EIGEN_HOME https://github.com/eigenteam/eigen-git-mirror/archive)
set(EIGEN_DOWNLOAD_URL ${EIGEN_HOME}/${EIGEN_VERSION}.tar.gz)
unset(EIGEN_HOME)

# Download and build the Eigen library and add its properties to the third party arguments.
set(EIGEN_ROOT ${THIRDPARTY_INSTALL_PATH} CACHE INTERNAL "")
ExternalProject_Add(eigen
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && curl -k -L ${EIGEN_DOWNLOAD_URL} -o eigen.tar.gz && tar -xzvf eigen.tar.gz && mv eigen-git-mirror-3.3.7 eigen && mv eigen ./Install
    URL ${EIGEN_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EIGEN_ROOT}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)

list(APPEND GLOBAL_THIRDPARTY_LIB_ARGS "-DEIGEN_ROOT:PATH=${EIGEN_ROOT}")
unset(EIGEN_DOWNLOAD_URL)
unset(EIGEN_ROOT)
