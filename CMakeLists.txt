project (crange) # VERSION 2.0.0)

cmake_minimum_required (VERSION 2.8)

if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif()

include(FeatureSummary)

set(${PROJECT_NAME}_VERSION_MAJOR 2)
set(${PROJECT_NAME}_VERSION_MINOR 0)
set(${PROJECT_NAME}_VERSION_PATCH 0)

# The version number.

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

add_library(crange_shared SHARED crange.cpp dictionary.c iniparser.c)
SET_TARGET_PROPERTIES(crange_shared PROPERTIES OUTPUT_NAME "CRange")
# add the executable
add_executable(crange main.cpp)
target_link_libraries(crange crange_shared)
add_dependencies(crange doc)

# Installation
FILE(GLOB headers "*.h")
FILE(GLOB readme "*.rst")
FILE(GLOB ini "*.ini")
install(FILES ${headers} DESTINATION include)
install(TARGETS crange crange_shared RUNTIME DESTINATION bin LIBRARY DESTINATION lib)
install(FILES ${readme} ${ini} tasks.txt DESTINATION share/doc/${PROJECT_NAME})
install(DIRECTORY html DESTINATION share/doc/${PROJECT_NAME})
