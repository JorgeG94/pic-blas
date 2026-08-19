if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "RelWithDebInfo"
      CACHE STRING "Build type to be used." FORCE)
  message(
    STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified.")
endif()

# specify module installation directory
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/modules/"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")

# Compiler-specific configurations
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(standard_compile_flags "-ffree-line-length-none -fbacktrace")
  set(warning_flags
      "-Wall -Wextra -Wpedantic -Wsurprising -Waliasing -Wampersand -Warray-bounds -Wcharacter-truncation -Wconversion -Wline-truncation -Wintrinsics-std -Wno-tabs -Wunderflow -Wunused-parameter -Wintrinsic-shadow -Wno-align-commons -Wsurprising"
  )
  set(fpe "-ffpe-trap=invalid,zero,overflow")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # -auto is required for thread safety: Intel's default -auto-scalar keeps
  # local scalars of INTEGER/REAL/COMPLEX/LOGICAL on the stack but leaves
  # CHARACTER locals in static storage. The transpose/uplo/side/jobz flags in
  # the BLAS and LAPACK wrappers are character(len=1) locals, so without -auto
  # they are shared by every thread and race whenever a caller runs these inside
  # an OpenMP region -- MKL then sees a TRANSA that disagrees with the LDA
  # computed from it and reports "Condition -1/-8 detected in function DGEMM".
  set(standard_compile_flags "-traceback -auto")
  set(fpe "-check all -fpe0")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  # see the -auto note above
  set(standard_compile_flags "-axAVX2 -traceback -auto")
  set(fpe "-check all -fpe0")
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NVHPC")
  set(standard_compile_flags "-traceback")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang")
  set(standard_compile_flags "")
endif()

# Customize compiler flags
set(CMAKE_Fortran_FLAGS_DEBUG
    "${CMAKE_Fortran_FLAGS_DEBUG} ${fpe}"
    PARENT_SCOPE)
set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} ${standard_compile_flags}"
    PARENT_SCOPE)

set(CMAKE_Fortran_FLAGS_COVERAGE
    "-O0 -g -coverage -fprofile-arcs -ftest-coverage")
if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
  set(CMAKE_Fortran_FLAGS
      "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_COVERAGE}"
      PARENT_SCOPE)
endif()
