cmake_minimum_required(VERSION 3.20)
project(PolyFit)


# c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# options
option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)


# headers
include_directories(PolyFit include)

# PolyFit library
add_library(
  PolyFit
  src/polyfit.cpp
)
install(TARGETS PolyFit)

# demo executable
add_executable(
  demo
  src/demo.cpp
)
target_link_libraries(
  demo
  PUBLIC
  PolyFit
)
