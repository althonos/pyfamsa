#.rst:
# FindPOPCNT
# --------
#
# Finds POPCNT support
#
# This module can be used to detect POPCNT support in a C compiler.  If
# the compiler supports POPCNT, the flags required to compile with
# POPCNT support are returned in variables for the different languages.
# The variables may be empty if the compiler does not need a special
# flag to support POPCNT.
#
# The following variables are set:
#
# ::
#
#    POPCNT_C_FLAGS - flags to add to the C compiler for POPCNT support
#    POPCNT_FOUND - true if POPCNT is detected
#
#=============================================================================

set(_POPCNT_REQUIRED_VARS)
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ${POPCNT_FIND_QUIETLY})

# sample POPCNT source code to test
set(POPCNT_C_TEST_SOURCE
"
#include <immintrin.h>
int foo(long long int a)
{
    return _popcnt64(a);
}
int main(void) { return (int)foo(0); }
")

# if these are set then do not try to find them again,
# by avoiding any try_compiles for the flags
if((DEFINED POPCNT_C_FLAGS) OR (DEFINED HAVE_POPCNT))
else()
  if(WIN32)
    # MSVC can compile AVX intrinsics without the arch flag, however it
    # will detect that AVX code is found and "consider using /arch:AVX".
    set(POPCNT_C_FLAG_CANDIDATES
      "/arch:SSE4.2"
    )
  else()
    set(POPCNT_C_FLAG_CANDIDATES
      #clang,GNU,Intel
      "-mpopcnt"
    )
  endif()

  include(CheckCSourceCompiles)

  foreach(FLAG IN LISTS POPCNT_C_FLAG_CANDIDATES)
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(HAVE_POPCNT CACHE)
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Try POPCNT C flag = [${FLAG}]")
    endif()
    check_c_source_compiles("${POPCNT_C_TEST_SOURCE}" HAVE_POPCNT)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(HAVE_POPCNT)
      set(POPCNT_C_FLAGS_INTERNAL "${FLAG}")
      break()
    endif()
  endforeach()

  unset(POPCNT_C_FLAG_CANDIDATES)
  
  set(POPCNT_C_FLAGS "${POPCNT_C_FLAGS_INTERNAL}"
    CACHE STRING "C compiler flags for POPCNT intrinsics")
endif()

list(APPEND _POPCNT_REQUIRED_VARS POPCNT_C_FLAGS)

set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if(_POPCNT_REQUIRED_VARS)
  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(POPCNT
                                    REQUIRED_VARS ${_POPCNT_REQUIRED_VARS})

  mark_as_advanced(${_POPCNT_REQUIRED_VARS})

  unset(_POPCNT_REQUIRED_VARS)
else()
  message(SEND_ERROR "FindPOPCNT requires C or CXX language to be enabled")
endif()
