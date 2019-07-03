if(APPLE)
file(GLOB dirs /Library/gurobi*)
foreach(d in ${dirs})
	string(REGEX MATCH "[0-9]+" GUROBI_VERSION "${d}")
endforeach(d)
elseif(UNIX)
file(GLOB dirs /opt/applications/gurobi/gurobi*)
foreach(d in ${dirs})
	string(REGEX MATCH "[0-9]+" GUROBI_VERSION "${d}")
endforeach(d)
endif()

message("Gurobi version ${GUROBI_VERSION}")
if(APPLE)
string(CONCAT GUROBI_DIR /Library/gurobi;${GUROBI_VERSION};/mac64)
elseif(UNIX)
string(CONCAT GUROBI_DIR /opt/applications/gurobi/gurobi;${GUROBI_VERSION};/linux64)
endif()
message("Looking for Gurobi in ${GUROBI_DIR}")

string(SUBSTRING ${GUROBI_VERSION} 0 2 GUROBI_VERSION_SHORT)

message("Gurobi version short ${GUROBI_VERSION_SHORT}")

find_path(GUROBI_INCLUDE_DIR gurobi_c++.h HINTS "${GUROBI_DIR}/include")
if(APPLE)
 find_library(GUROBI_LIBRARY libgurobi${GUROBI_VERSION_SHORT}.dylib HINTS ${GUROBI_DIR}/lib)
elseif(UNIX)
 find_library(GUROBI_LIBRARY libgurobi${GUROBI_VERSION_SHORT}.so HINTS ${GUROBI_DIR}/lib)
endif()
find_library(GUROBI_CPP_LIBRARY libgurobi_c++.a HINTS ${GUROBI_DIR}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GUROBI DEFAULT_MSG GUROBI_LIBRARY GUROBI_CPP_LIBRARY GUROBI_INCLUDE_DIR)

if(GUROBI_FOUND)
    set(GRB_LICENSE_FILE "~/gurobi.research.lic")
    set(GUROBI_INCLUDE_DIRS ${GUROBI_INCLUDE_DIR})
    set(GUROBI_LIBRARIES ${GUROBI_CPP_LIBRARY} ${GUROBI_LIBRARY})
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        set(GUROBI_LIBRARIES "${GUROBI_LIBRARIES};m;pthread")
    endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(GUROBI_FOUND)

mark_as_advanced(GUROBI_LIBRARY GUROBI_CPP_LIBRARY GUROBI_INCLUDE_DIR)
