cmake_minimum_required(VERSION 3.16)
project(swpc_psv_cpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(SWPC_OPENMP "Enable OpenMP threading of the FDM kernels" ON)
# Strict IEEE mode (no fast-math, no FMA contraction) is the default so that
# differential validation against the Fortran reference is reproducible.
option(SWPC_FAST_MATH "Allow -ffast-math (breaks bitwise reproducibility)" OFF)

file(GLOB SWPC_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_executable(swpc_psv_cpp ${SWPC_SOURCES})
target_compile_options(swpc_psv_cpp PRIVATE -Wall -Wextra -ffp-contract=off)
if(SWPC_FAST_MATH)
  target_compile_options(swpc_psv_cpp PRIVATE -ffast-math)
endif()
if(SWPC_OPENMP)
  find_package(OpenMP REQUIRED)
  target_link_libraries(swpc_psv_cpp PRIVATE OpenMP::OpenMP_CXX)
endif()

enable_testing()
