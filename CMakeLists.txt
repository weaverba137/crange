project (Tutorial)

cmake_minimum_required (VERSION 3.0)

include(FeatureSummary)

# The version number.
set (Tutorial_VERSION_MAJOR 1)
set (Tutorial_VERSION_MINOR 0)
set (Tutorial_VERSION_PATCH 1)
set (Tutorial_VERSION "${Tutorial_VERSION_MAJOR}.${Tutorial_VERSION_MINOR}.${Tutorial_VERSION_PATCH}")
set (Tutorial_VERSION_STRING "${Tutorial_VERSION}")

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
        "${PROJECT_SOURCE_DIR}/TutorialConfig.h.in"
        "${PROJECT_BINARY_DIR}/TutorialConfig.h"
    )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
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
        "${PROJECT_SOURCE_DIR}/Doxyfile.cmake.in"
        "${PROJECT_BINARY_DIR}/Doxyfile" @ONLY
    )
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

add_library(TutorialFunctions SHARED tutorial.cpp dictionary.c iniparser.c)
# add the executable
add_executable(Tutorial main.cpp)
target_link_libraries(Tutorial TutorialFunctions)
