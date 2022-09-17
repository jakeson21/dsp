cmake_minimum_required(VERSION 3.1)
project(dsp)

set(PROJECT_VERSION 1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#######################
add_library(dsp
    filters/iir_bandpass.cpp
    filters/iir_bandstop.cpp
)

target_include_directories(dsp PUBLIC
    filters)

#######################
add_executable(dsp-test 
    test/dsp-test.cpp)

target_include_directories(dsp-test PUBLIC
    test
    filters)

target_link_libraries(dsp-test PUBLIC
    dsp)

# find_package(Matplot++ QUIET)
# if(NOT Matplot++_FOUND)
#     message(STATUS "Matplot++ package not found, fetching it now...")
#     include(FetchContent)
#     FetchContent_Declare(matplotplusplus
#         GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
#         GIT_TAG origin/master) # or whatever tag you want
#     FetchContent_GetProperties(matplotplusplus)
#     if(NOT matplotplusplus_POPULATED)
#         FetchContent_Populate(matplotplusplus)
#         add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
#     endif()

#     message(STATUS "Matplot++ has been fetched...")

#     target_link_libraries(dsp-test PUBLIC matplot)
# endif()

# For SigPack
find_package(Armadillo 7.6 REQUIRED)
message(STATUS "ARMADILLO_LIBRARIES: ${ARMADILLO_LIBRARIES}")
message(STATUS "ARMADILLO_INCLUDE_DIRS: ${ARMADILLO_INCLUDE_DIRS}")
target_link_libraries(dsp-test PUBLIC ${ARMADILLO_LIBRARIES})
target_include_directories(dsp-test PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/extern/SigPack/sigpack/
    ${ARMADILLO_INCLUDE_DIRS}
)
