cmake_minimum_required( VERSION 2.8 )

#################################################
## Compilateur
#################################################

set (CMAKE_CXX_STANDARD 11)

## C++ compiler options
set( CMAKE_CXX_FLAGS "-Wall -Wextra" )
set( CMAKE_CXX_FLAGS_DEBUG "-g -O0" )
set( CMAKE_CXX_FLAGS_RELEASE "-O3" )

#################################################
## Project
#################################################

## Sources and headers
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

## Find all source files to be compiled into object files
file( GLOB_RECURSE SOURCES *.cc *.cpp )

## Source file containing the "main" function
set( MAIN_SOURCES main.cpp )

## Specify a name for the generated executable file
set( MAIN_EXE_NAME Main )

## 1. Compile...
add_executable( ${MAIN_EXE_NAME} ${MAIN_SOURCES} ${SOURCES})

#################################################
## Depedancies
#################################################

# Visp
find_package(VISP REQUIRED)
include_directories(${VISP_INCLUDE_DIRS})

target_link_libraries( ${MAIN_EXE_NAME} ${VISP_LIBRARIES})