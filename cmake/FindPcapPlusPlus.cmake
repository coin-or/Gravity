# FindPcapPlusPlus.cmake
#
# Finds the PcapPlusPlus library.
#
# This will define the following variables
#
#    PcapPlusPlus_FOUND
#    PcapPlusPlus_INCLUDE_DIRS
#    PcapPlusPlus_LIBRARIES
#    PcapPlusPlus_VERSION
#
# and the following imported targets
#
#    PcapPlusPlus::PcapPlusPlus
#

if (PC_PcapPlusPlus_INCLUDEDIR AND PC_PcapPlusPlus_LIBDIR)
    set(PcapPlusPlus_FIND_QUIETLY TRUE)
endif ()

find_package(PkgConfig REQUIRED)
pkg_check_modules(PC_PcapPlusPlus REQUIRED PcapPlusPlus)

set(PcapPlusPlus_VERSION ${PC_PcapPlusPlus_VERSION})

mark_as_advanced(PcapPlusPlus_INCLUDE_DIR PcapPlusPlus_LIBRARY)

foreach (LIB_NAME ${PC_PcapPlusPlus_LIBRARIES})
    find_library(${LIB_NAME}_PATH ${LIB_NAME} HINTS ${PC_PcapPlusPlus_LIBDIR})
    if (${LIB_NAME}_PATH)
        list(APPEND PcapPlusPlus_LIBS ${${LIB_NAME}_PATH})
    endif ()
endforeach ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PcapPlusPlus
        REQUIRED_VARS PC_PcapPlusPlus_INCLUDEDIR PC_PcapPlusPlus_LIBDIR
        VERSION_VAR PcapPlusPlus_VERSION
        )

if (PcapPlusPlus_FOUND)
    set(PcapPlusPlus_INCLUDE_DIRS ${PC_PcapPlusPlus_INCLUDEDIR})
    set(PcapPlusPlus_LIBRARIES ${PcapPlusPlus_LIBS})
endif ()

if (PcapPlusPlus_FOUND AND NOT TARGET PcapPlusPlus::PcapPlusPlus)
    add_library(PcapPlusPlus::PcapPlusPlus INTERFACE IMPORTED)
    set_target_properties(PcapPlusPlus::PcapPlusPlus PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${PcapPlusPlus_INCLUDE_DIRS}"
            INTERFACE_LINK_LIBRARIES "${PcapPlusPlus_LIBRARIES}"
            INTERFACE_COMPILE_FLAGS "${PC_PcapPlusPlus_CFLAGS}"
            )
endif ()
