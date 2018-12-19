cmake_minimum_required(VERSION 3.5)
cmake_policy(SET CMP0074 NEW)

project(GeoCARBSIF)
set(PROJECT_DESCRIPTION "GeoCARB SIF Retrieval Algorithm")
enable_language(Fortran)

# This contains extra CMAKE modules to find NetCDF libraries
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# This here lets us use the current GIT revision as part of
# the make process (TODO!!)

execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

execute_process(
  COMMAND git rev-list --count HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_REV_NO
  OUTPUT_STRIP_TRAILING_WHITESPACE
  )

add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
add_definitions("-DGIT_BRANCH=${GIT_BRANCH}")
add_definitions("-DGIT_REV_NO=${GIT_REV_NO}")

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parameters/Version.f90.in
  ${CMAKE_CURRENT_BINARY_DIR}/generated/Version.f90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(dialect "-ffree-form -std=f2018 -fimplicit-none")
  # GCC writes a lot of bogus warnings when dealing with allocatable arrays, so
  # have to turn those warnings off.
  set(warnings "-Wall -Wextra -Wpedantic -Wno-uninitialized")
  set(opt "-O3")
  set(bounds "-fbounds-check")
else()
  # only support Gfortran for now
endif()

if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 8.0)
    message(FATAL_ERROR "gfortran version must be at least 8.0!")
else()
    message(STATUS "gfortran version is sufficcient!")
endif()

set(CMAKE_Fortran_FLAGS "${dialect}")
set(CMAKE_Fortran_FLAGS_DEBUG "-g ${bounds} ${warnings}")
set(CMAKE_Fortran_FLAGS_RELEASE "${opt}")

# Here are the standard modules we need in order to compile
# properly. Most machines should have these by now.
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS Fortran HL)
find_package(NETCDF REQUIRED COMPONENTS Fortran)

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
# all files in the third-party directory. Also, most third party stuff has not
# been written to the latest standard, or best practices..
set_source_files_properties(
  ${third-party}
  PROPERTIES
  COMPILE_FLAGS "-w"
)

# Define the executable
add_executable(GeoCARBSIF ${sources} ${third-party} ${generated})

message(STATUS "HDF5-HL Libraries at -- " ${HDF5_Fortran_HL_LIBRARIES})
message(STATUS "HDF5 Libraries at -- " ${HDF5_Fortran_LIBRARIES})

# HDF5 needs to be included at compile time already
target_include_directories(GeoCARBSIF PUBLIC ${HDF5_Fortran_INCLUDE_DIRS} ${HDF5_Fortran_HL_INCLUDE_DIRS})
target_link_libraries(GeoCARBSIF ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES} ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES})
