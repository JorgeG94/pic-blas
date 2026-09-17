cmake_minimum_required(VERSION 3.22)
cmake_policy(SET CMP0144 NEW)

# A minimal consumer of an *installed* pic-blas, built by CI only.
#
# The FetchContent consumer next door covers the subproject path. This one
# covers the other one: find_package against an install tree. Nothing else in CI
# runs `cmake --install` at all, so pic-blasConfig.cmake was never once
# exercised -- it could name dependencies that do not resolve, or omit ones it
# needs, and every job would still be green.
project(
  "pic-blas-installed-consumer"
  LANGUAGES Fortran C
  VERSION 0.1)

enable_testing()

# No CMAKE_MODULE_PATH, deliberately: a real consumer does not have pic-blas's
# cmake/modules/Findpic.cmake, so this resolves pic the way they would -- via
# the installed picConfig.cmake that our own config has to ask for.
find_package(pic-blas REQUIRED)

if(NOT TARGET pic-blas::pic-blas)
  message(
    FATAL_ERROR "find_package(pic-blas) did not define pic-blas::pic-blas")
endif()

# Same source as the FetchContent consumer: one program, both consumption paths.
add_executable(installed_consumer
               "${CMAKE_CURRENT_LIST_DIR}/../consumer/main.f90")
target_link_libraries(installed_consumer PRIVATE pic-blas::pic-blas)

add_test(NAME installed_consumer_runs COMMAND installed_consumer)
