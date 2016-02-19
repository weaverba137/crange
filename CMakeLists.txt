project (CRange)

cmake_minimum_required (VERSION 2.8)

include(FeatureSummary)

# The version number.
set (crange_VERSION_MAJOR 2)
set (crange_VERSION_MINOR 0)
set (crange_VERSION_PATCH 0)
set (crange_VERSION "${crange_VERSION_MAJOR}.${crange_VERSION_MINOR}.${crange_VERSION_PATCH}")
set (crange_VERSION_STRING "${crange_VERSION}")

set (CMAKE_CXX_FLAGS "-std=c++11 -Wall")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
        "${PROJECT_SOURCE_DIR}/config.h.in"
        "${PROJECT_BINARY_DIR}/config.h"
    )

# add the binary tree to the search path for include files
# so that we will find config.h
include_directories("${PROJECT_BINARY_DIR}")

# Doxygen support
find_package(Doxygen)
set_package_properties(Doxygen PROPERTIES
    TYPE OPTIONAL
    DESCRIPTION "Doxygen documentation generator"
    PURPOSE "Needed for generating API documentation (make doc)"
)
if(DOXYGEN_FOUND)
    configure_file(
        "${PROJECT_SOURCE_DIR}/doxygen.cfg.in"
        "${PROJECT_BINARY_DIR}/doxygen.cfg" @ONLY
    )
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doxygen.cfg
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

add_library(CRange SHARED crange.cpp dictionary.c iniparser.c)
# add the executable
add_executable(crange main.cpp)
target_link_libraries(crange CRange)
