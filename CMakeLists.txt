cmake_minimum_required(VERSION 3.1)
project(dsp)

set(PROJECT_VERSION 1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#######################
add_library(dsp
    filters/iir_bandpass.cpp)

target_include_directories(dsp PUBLIC
    filters)

#######################
add_executable(dsp-test 
    test/dsp-test.cpp)

target_include_directories(dsp-test PUBLIC
    test
    filters)

target_link_libraries(dsp-test 
    dsp)
