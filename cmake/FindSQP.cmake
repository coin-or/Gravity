set(RestartSQP_ROOT_DIR "$ENV{RestartSQP_ROOT_DIR}" CACHE PATH "RestartSQP root directory.")
message("Looking for RestartSQP in ${RestartSQP_ROOT_DIR}")
set(QPOASES_ROOT_DIR "$ENV{QPOASES_ROOT_DIR}" CACHE PATH "QPOASES root directory.")
if(QPOASES_ROOT_DIR)
 message("Looking for QPOASES in ${QPOASES_ROOT_DIR}")
else(QPOASES_ROOT_DIR)
 message("QPOASES_ROOT_DIR not provided.")
endif(QPOASES_ROOT_DIR)

set(QORE_ROOT_DIR "$ENV{QORE_ROOT_DIR}" CACHE PATH "QORE root directory.")
message("Looking for QORE in ${QORE_ROOT_DIR}")

find_path(RestartSQP_INCLUDE_DIR
	NAMES LazySqpSolver.hpp 
	HINTS /usr/local/include
	HINTS ${RestartSQP_ROOT_DIR}/include/RestartSQP
	HINTS ${RestartSQP_ROOT_DIR}/../include/restartsqp
	HINTS ${RestartSQP_ROOT_DIR}/include/restartsqp
)

find_library(QPOASES_LIBRARY
  qpOASES
  HINTS ${QPOASES_ROOT_DIR}/build/libs
  HINTS ${QPOASES_ROOT_DIR}/bin
  HINTS ${QPOASES_ROOT_DIR}/lib
  HINTS ${QPOASES_ROOT_DIR}/libs
  HINTS ${QPOASES_ROOT_DIR}/lib/Release
)


find_library(QORE_LIBRARY
        qore
        HINTS /usr/local/lib
        HINTS third_party/QORE/lib
        HINTS ${QORE_ROOT_DIR}/test_build/lib
        HINTS ${QORE_ROOT_DIR}/lib
        HINTS ${QORE_ROOT_DIR}/build/lib/Release
)

find_library(RestartSQP_LIBRARY 
	librestartsqp.a
	HINTS /usr/local/lib
	HINTS ${RestartSQP_ROOT_DIR}/lib
	HINTS ${RestartSQP_ROOT_DIR}-build/lib
	HINTS ${RestartSQP_ROOT_DIR}/libs
	HINTS ${RestartSQP_ROOT_DIR}/lib/Release
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(RestartSQP DEFAULT_MSG RestartSQP_LIBRARY RestartSQP_INCLUDE_DIR)

if(RestartSQP_FOUND)
	message("—- Found RestartSQP header files under ${RestartSQP_INCLUDE_DIR}")
	message("—- Found RestartSQP library under ${RestartSQP_LIBRARY}")
    set(SQP_INCLUDE_DIRS ${RestartSQP_INCLUDE_DIR})
    set(SQP_LIBRARIES ${RestartSQP_LIBRARY} ${QORE_LIBRARY} ${QPOASES_LIBRARY} -lblas -llapack)
endif(RestartSQP_FOUND)

mark_as_advanced(RestartSQP_LIBRARY RestartSQP_INCLUDE_DIR)
