cmake_minimum_required (VERSION 2.8.11)

project (templates)

add_executable (square square.cc)
add_executable (functiontemplate functiontemplate.cc)
add_executable (classes classes.cc)
add_executable (explicit explicit.cc)

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-Wall")
endif()

