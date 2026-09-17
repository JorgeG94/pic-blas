cmake_minimum_required(VERSION 3.22)
cmake_policy(SET CMP0144 NEW)

# A minimal downstream consumer, built by CI only.
#
# pic-blas's own jobs build it standalone, where CMAKE_BINARY_DIR and
# PROJECT_BINARY_DIR are the same directory -- so a target that leaks the
# consumer's build tree is invisible there. This project consumes pic-blas the
# way a real downstream project does, via FetchContent, where the two differ. It
# exists to fail if that ever regresses.
project(
  "pic-blas-consumer"
  LANGUAGES Fortran C
  VERSION 0.1)

enable_testing()

include(FetchContent)

# A consumer wants the library, not the dependency chain's test suites.
# ENABLE_TESTING is pic-blas's own switch; PIC_ENABLE_TESTING is pic's, and
# without it pic's 40-odd tests get built and run inside this job.
set(ENABLE_TESTING
    OFF
    CACHE BOOL "" FORCE)
set(PIC_ENABLE_TESTING
    OFF
    CACHE BOOL "" FORCE)

FetchContent_Declare(pic-blas SOURCE_DIR "${PIC_BLAS_SOURCE_DIR}")
FetchContent_MakeAvailable(pic-blas)

# Mirror the fallback in a consumer's fetch helper: if pic-blas does not provide
# the namespaced target itself, the consumer invents an INTERFACE IMPORTED one.
# CMake validates INTERFACE_INCLUDE_DIRECTORIES on imported targets, so this is
# the path that turns a wrong module directory into a hard configure error.
if(NOT TARGET pic-blas::pic-blas)
  message(WARNING "pic-blas did not provide pic-blas::pic-blas; inventing one")
  add_library(pic-blas::pic-blas INTERFACE IMPORTED)
  target_link_libraries(pic-blas::pic-blas INTERFACE pic-blas)
endif()

add_executable(consumer main.f90)
target_link_libraries(consumer PRIVATE pic-blas::pic-blas)

add_test(NAME consumer_runs COMMAND consumer)

# pic-blas must build into its own tree, not ours.
if(EXISTS "${CMAKE_BINARY_DIR}/modules")
  message(
    FATAL_ERROR
      "pic-blas wrote its module directory into the consumer's build tree "
      "(${CMAKE_BINARY_DIR}/modules). It must use its own PROJECT_BINARY_DIR.")
endif()
