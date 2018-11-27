cmake_minimum_required(VERSION 3.5)

project(GeoCARBSIF)
set(PROJECT_DESCRIPTION "GeoCARB SIF Retrieval Algorithm")
enable_language(Fortran)

# This contains extra CMAKE modules to find NetCDF libraries
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# This here lets us use the current GIT revision as part of
# the make process (TODO!!)

execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
add_definitions("-DGIT_BRANCH=${GIT_BRANCH}")

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parameters/Version.F90.in
  ${CMAKE_CURRENT_BINARY_DIR}/generated/Version.F90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-ffree-form -std=f2008 -fimplicit-none")
  # GCC writes a lot of bogus warnings when dealing with allocatable arrays, so
  # have to turn those warnings off.
  set(warnings "-Wall -Wextra -Wpedantic -Wno-uninitialized")
  set(opt "-O3")
  set(bounds "-fbounds-check")
endif()

set(CMAKE_Fortran_FLAGS "${dialect}")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 ${bounds} ${warnings}")
set(CMAKE_Fortran_FLAGS_RELEASE "${opt}")

# Here are the standard modules we need in order to compile
# properly. Most machines should have these by now.
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(HDF5 REQUIRED)
find_package(NETCDF REQUIRED)

# This is where we put compiled modules
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Add source files to the list
file(GLOB_RECURSE sources src/*.?90)
# and the third party libraries
file(GLOB_RECURSE third-party third-party/*.?90)
# Add generated files to the list
file(GLOB_RECURSE generated ${CMAKE_CURRENT_BINARY_DIR}/generated/*.?90)

# Now we want warnings, but not for third-party libraries,
# we sort of trust those to be working. So just disable warnings for
# all files in the third-party directory.
set_source_files_properties(
  ${third-party}
  PROPERTIES
  COMPILE_FLAGS "-w"
)

# Define the executable
add_executable(GeoCARBSIF ${sources} ${third-party} ${generated})
# Link the BLAS and LAPACK libraries here
target_link_libraries(GeoCARBSIF lapack blas)
