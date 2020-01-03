project(eigen)

cmake_minimum_required(VERSION 2.6)

find_package(VISP REQUIRED)

find_package( OpenMP REQUIRED )
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(cpp_
  main.cpp )

foreach(cpp ${cpp_})
  visp_add_target(${cpp})
endforeach()

target_link_libraries( main OpenMP::OpenMP_CXX)
