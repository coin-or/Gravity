# Create download URL derived from version number.
set(OPT_PARSER_HOME https://github.com/Lcressot/cpp_option_parser/archive)
set(OPT_PARSER_DOWNLOAD_URL ${OPT_PARSER_HOME}/${OPT_PARSER_VERSION}.tar.gz)
unset(OPT_PARSER_HOME)

# Download and build the Eigen library and add its properties to the third party arguments.
set(OPT_PARSER_ROOT ${THIRDPARTY_INSTALL_PATH} CACHE INTERNAL "")
ExternalProject_Add(opt_parser
    DOWNLOAD_DIR ${THIRDPARTY_INSTALL_PATH}
    DOWNLOAD_COMMAND export HTTPS_PROXY=$ENV{HTTPS_PROXY} && curl -k -L ${OPT_PARSER_DOWNLOAD_URL} -o cpp_option_parser.tar.gz && tar -xzvf cpp_option_parser.tar.gz && mv cpp_option_parser-1.0 opt_parser && rm -fr ./Source/opt_parser && mv ./opt_parser ./Source
    URL ${OPT_PARSER_DOWNLOAD_URL}
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${OPT_PARSER_ROOT}
    INSTALL_COMMAND ""
    UPDATE_COMMAND ""
)

unset(OPT_PARSER_DOWNLOAD_URL)
unset(OPT_PARSER_ROOT)
