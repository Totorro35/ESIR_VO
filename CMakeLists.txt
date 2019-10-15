cmake_minimum_required( VERSION 2.8 )

## Load CMAKE configuration from environment variables
set( CMAKE_MODULE_PATH $ENV{CMAKE_MODULE_PATH} )
set( CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH} )

#project name need to be in 1 word
project(TP_VO)

add_subdirectory(${PROJECT_SOURCE_DIR}/src bin)
add_subdirectory(${PROJECT_SOURCE_DIR}/docs docs)
add_subdirectory(${PROJECT_SOURCE_DIR}/tp tp)