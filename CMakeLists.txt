cmake_minimum_required(VERSION 3.22)
cmake_policy(SET CMP0144 NEW)

project(
  "pic-blas"
  LANGUAGES Fortran C
  VERSION 0.3.3
  DESCRIPTION "ADD YOUR DESCRIPTION HERE")

set(project_name pic-blas)
set(main_lib ${project_name})
set(all_targets ${main_lib})

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(GNUInstallDirs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

# see cmake/ subdirectory for this
include(demo_status)

option(ENABLE_TESTING "Enable testing for ${project_name}" ON)

add_subdirectory(cmake)

set(libraries_to_link)
if(NOT TARGET test-drive::test-drive AND ENABLE_TESTING)
  enable_testing()
  find_package("test-drive" REQUIRED)
endif()

if(NOT TARGET pic::pic)
  find_package("pic" REQUIRED)
endif()

add_library(${main_lib} STATIC)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# pic is PUBLIC, not PRIVATE: pic_blas_interfaces and pic_lapack_interfaces both
# `use pic_types`, so pic is part of this library's public interface -- a
# consumer compiling against our .mod files needs pic's too. PRIVATE hid pic's
# module directory from everything downstream, including our own test target,
# and only worked while pic happened to write its .mod files into the top-level
# build tree; pic 0.8.2 moved them into its own, which broke that accident.
target_link_libraries(${main_lib} PUBLIC pic::pic)
target_link_libraries(${main_lib} PRIVATE BLAS::BLAS LAPACK::LAPACK)

# Enable position-independent code for shared-library linking support
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(src)

if(ENABLE_TESTING)
  add_subdirectory(test)
endif()

set_target_properties(${main_lib} PROPERTIES Fortran_MODULE_DIRECTORY
                                             ${CMAKE_BINARY_DIR}/modules)

# Carry our own module directory to anyone who links us, so no consumer has to
# hardcode where the .mod files land -- the same spelling that just broke.
target_include_directories(
  ${main_lib} INTERFACE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>)

foreach(tgt ${all_targets})
  message(STATUS "${tgt}")
  install(
    TARGETS ${tgt}
    EXPORT ${project_name}Targets
    ARCHIVE DESTINATION lib
    INCLUDES
    DESTINATION include/${project_name})
endforeach()

install(
  DIRECTORY ${CMAKE_BINARY_DIR}/modules/
  DESTINATION include/${project_name}
  FILES_MATCHING
  PATTERN "*.mod")

install(
  EXPORT ${project_name}Targets
  FILE ${project_name}Targets.cmake
  NAMESPACE ${project_name}::
  DESTINATION lib/cmake/${project_name})

include(CMakePackageConfigHelpers)

# RENAME YOUR sampleConfig.cmake.in to match ${project_name}Config.cmake
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/${project_name}Config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/${project_name}Config.cmake"
  INSTALL_DESTINATION lib/cmake/${project_name})

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${project_name}ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${project_name}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${project_name}ConfigVersion.cmake"
        DESTINATION lib/cmake/${project_name})

install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/${project_name}Config.cmake"
  DESTINATION lib/cmake/${project_name}
  RENAME ${project_name}-config.cmake)

print_summary()

if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
  add_custom_target(
    coverage
    COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
    COMMAND ${CMAKE_COMMAND} -E remove -f coverage.info
    COMMAND lcov --directory . --capture --output-file coverage.info
    COMMAND lcov --ignore-errors unused --remove coverage.info '/usr/*'
            '*/build/*' --output-file coverage_filtered.info
    COMMAND genhtml coverage_filtered.info --output-directory coverage_report
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating code coverage report")
endif()
