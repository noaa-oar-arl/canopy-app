# Top-level CMakeLists.txt for canopy-app

cmake_minimum_required(VERSION 3.15)
project(canopy-app LANGUAGES Fortran)

# Optionally set a default build type if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Add the src directory (where the main CMakeLists.txt and sources are)
add_subdirectory(src)
