project(tp_av)

cmake_minimum_required(VERSION 2.6)

find_package(VISP REQUIRED)

set(cpp_
  av-v0.cpp
  eigenfaces.cpp)

foreach(cpp ${cpp_})
  visp_add_target(${cpp})
endforeach()
