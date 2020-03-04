cmake_minimum_required(VERSION 3.5)
cmake_policy(SET CMP0074 NEW)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")

message(STATUS "CMake module path is: ${CMAKE_MODULE_PATH}")

project(GASBAG)

set(PROJECT_DESCRIPTION "Generic Algorithm for Single Band Acquisition of Gases - GASBAG")
enable_language(Fortran C)


# This section adds the git branch, commit hash and "revision number"
# as compile-time variables. So running ./GASBAG will print out those
# values.

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


# We have currently these build types. DEBUG will
# add the "DEBUG" preprocessor varible.

string( TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE )

if (BUILD_TYPE STREQUAL "debug")
  message(STATUS "Compiling in debug mode.")
  add_definitions("-DDEBUG")
elseif (BUILD_TYPE STREQUAL "release")
  message(STATUS "Compiling in release mode.")
elseif (BUILD_TYPE STREQUAL "profile")
  message(STATUS "Compiling in profile mode.")
else()
  message(FATAL_ERROR "Build type ${BUILD_TYPE} not supported.
                       (DEBUG/RELEASE/PROFILE, case-insenstive)")
endif()

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/src/parameters/Version.f90.in
  ${CMAKE_CURRENT_BINARY_DIR}/generated/Version.f90
)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")

  if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 8.0)
      message(STATUS "gfortran version is sufficcient!")
  else()
      message(FATAL_ERROR "gfortran version must be 8.2! 9.2 produces problems!")
  endif()

  # MCMODEL LARGE MUST BE TURNED ON!
  # OTHERWISE YOU GET STRANGE ERRORS IN AN OPENMP ENVIRONMENT
  set(dialect "-cpp -ffree-form -std=f2018 -fimplicit-none -ffree-line-length-512")

  # Make OpenMP conditional
  if(USE_OPENMP)
    set(dialect "${dialect} -fopenmp")
    message(STATUS "Compiling with OpenMP")
  endif()

  # GCC writes a lot of bogus warnings when dealing with allocatable arrays, so
  # have to turn those warnings off.
  #set(warnings "-Wall -Wextra -Wpedantic -Wno-uninitialized")
  set(warnings "-Wall -Wextra -Wpedantic")
  set(opt "-O2 -mtune=native -march=native")
  set(debug "-Og -g3 -fcheck=all -fbacktrace -fno-omit-frame-pointer")
  set(profile "-Og -pg -fbacktrace -fprofile-generate --coverage")

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")

  set(dialect "-assume realloc_lhs -standard-semantics -std18 -debug all -warn all")
  if(USE_OPENMP)
    set(dialect "${dialect}")
    message(STATUS "Compiling with OpenMP")
  endif()

else()
  message(FATAL_ERROR "Unsupported compiler")
endif()



set(CMAKE_Fortran_FLAGS "${dialect}")
set(CMAKE_Fortran_FLAGS_DEBUG "${debug} ${warnings}")
set(CMAKE_Fortran_FLAGS_RELEASE "${opt}")
set(CMAKE_Fortran_FLAGS_PROFILE "${profile}")

set(CMAKE_CXX_FLAGS "-O2 -ffree-line-length-0")
#set(CMAKE_CXX_FLAGS_PROFILE "-Og -pg -fbacktrace -fprofile-arcs -ftest-coverage")

# Here are the standard modules we need in order to compile
# properly. Most machines should have these by now.
# NOTE OpenBLAS seems to be pretty quick. Outsourcing many operations
# to BLAS-equivalents might speed up the code even further.

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(NETCDF REQUIRED COMPONTENTS Fortran)
find_package(HDF5 REQUIRED COMPONENTS Fortran)
# This is where we put compiled modules
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)

# Add source files to the list
file(GLOB_RECURSE sources src/*.?90)
# and the third party libraries
file(GLOB_RECURSE third-party third-party/*.?90)
# Add generated files to the list
file(GLOB_RECURSE generated ${CMAKE_CURRENT_BINARY_DIR}/generated/*.?90)

# Add Greg's XRTM
# This is a little bit of a manual hackjob, and one can argue that
# using Greg's makefile setup might be nice. One would compile XRTM
# separately and just link to the resulting library. This however, lets
# you make changes to compiler flags more easily etc.

set(xrtm_sources
  ####
third-party/xrtm/src/version.c
third-party/xrtm/src/xrtm_adding_a.c
third-party/xrtm/src/xrtm_adding.c
third-party/xrtm/src/xrtm_brdf_a.c
third-party/xrtm/src/xrtm_brdf.c
third-party/xrtm/src/xrtm_derivs.c
third-party/xrtm/src/xrtm_doubling.c
third-party/xrtm/src/xrtm_doub_rts.c
third-party/xrtm/src/xrtm_eig_bvp_a.c
third-party/xrtm/src/xrtm_eig_bvp.c
third-party/xrtm/src/xrtm_eig_rts_a.c
third-party/xrtm/src/xrtm_eig_rts.c
third-party/xrtm/src/xrtm_eig_util_a.c
third-party/xrtm/src/xrtm_eig_util.c
third-party/xrtm/src/xrtm_external.c
third-party/xrtm/src/xrtm_four_stream.c
third-party/xrtm/src/xrtm_interface.c
third-party/xrtm/src/xrtm_matrix.c
third-party/xrtm/src/xrtm_mem_bvp.c
third-party/xrtm/src/xrtm_model_a.c
third-party/xrtm/src/xrtm_model.c
third-party/xrtm/src/xrtm_pade_rts.c
third-party/xrtm/src/xrtm_radiance_a.c
third-party/xrtm/src/xrtm_radiance.c
third-party/xrtm/src/xrtm_save_tree.c
third-party/xrtm/src/xrtm_scatter_a.c
third-party/xrtm/src/xrtm_scatter.c
third-party/xrtm/src/xrtm_sfi.c
third-party/xrtm/src/xrtm_single_a.c
third-party/xrtm/src/xrtm_single.c
third-party/xrtm/src/xrtm_six_stream.c
third-party/xrtm/src/xrtm_sos.c
third-party/xrtm/src/xrtm_source_a.c
third-party/xrtm/src/xrtm_source.c
third-party/xrtm/src/xrtm_stacks.c
third-party/xrtm/src/xrtm_support.c
third-party/xrtm/src/xrtm_two_os.c
third-party/xrtm/src/xrtm_two_stream.c
third-party/xrtm/src/xrtm_utility_a.c
third-party/xrtm/src/xrtm_utility.c
third-party/xrtm/src/xrtm_work.c
  ####
third-party/xrtm/misc/gindex_name_value.c
third-party/xrtm/misc/glist.c
third-party/xrtm/misc/gmath_mat2d3d.c
third-party/xrtm/misc/gmath_matrix.c
third-party/xrtm/misc/gmath_vector.c
third-party/xrtm/misc/gstruct.c
third-party/xrtm/misc/gutil_alloc.c
third-party/xrtm/misc/gutil_cmp.c
third-party/xrtm/misc/gutil_error.c
third-party/xrtm/misc/gutil_parse.c
third-party/xrtm/misc/gutil_sort.c
third-party/xrtm/misc/gutil_string.c
third-party/xrtm/misc/rtutil_math.c
third-party/xrtm/misc/rtutil_scat_fft.c
third-party/xrtm/misc/rtutil_scat_io.c
third-party/xrtm/misc/rtutil_size_dist.c
third-party/xrtm/misc/rtutil_support.c
third-party/xrtm/misc/zeroin.c
  ####
third-party/xrtm/interfaces/xrtm_int_f90_util.c
          )



# Third-party libraries are assumed to work fine, so compile with optimisation
# and no debugging options.
set_source_files_properties(
  ${third-party}
  PROPERTIES
  COMPILE_FLAGS "${opt}"
  COMPILE_FLAGS_DEBUG ""
)

set_source_files_properties(
  ${xrtm_sources}
  PROPERTIES
  COMPILE_FLAGS "-O2"
  COMPILE_FLAGS_DEBUG ""
  )

# Our code is definitely in the need of having warnings
# show up!
set_source_files_properties(
  ${sources}
  PROPERTIES
  COMPILE_FLAGS "${warnings}"
)


# Define the executable
add_executable(GASBAG ${sources} ${third-party} ${xrtm_sources} ${generated})

target_include_directories(GASBAG PRIVATE
  "third-party/xrtm"
  "third-party/xrtm/interfaces"
  "third-party/xrtm/utils"
  "third-party/xrtm/misc"
  "third-party/xrtm/src")

message(STATUS "LAPACK Libraries at -- " ${LAPACK_LIBRARIES})
message(STATUS "BLAS Libraries at -- " ${BLAS_LIBRARIES})
message(STATUS "HDF5 Libraries at -- " ${HDF5_Fortran_LIBRARIES})

# HDF5 needs to be included at compile time already
target_include_directories(GASBAG PUBLIC ${HDF5_Fortran_INCLUDE_DIRS}
#${GPERFTOOLS_INCLUDE_DIR}
)

target_link_libraries(GASBAG
${LAPACK_LIBRARIES}
${BLAS_LIBRARIES}
${HDF5_Fortran_LIBRARIES}
)
