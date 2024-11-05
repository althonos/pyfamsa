#.rst:
# FindAVX1
# --------
#
# Finds AVX1 support
#
# This module can be used to detect AVX1 support in a C compiler.  If
# the compiler supports AVX1, the flags required to compile with
# AVX1 support are returned in variables for the different languages.
# The variables may be empty if the compiler does not need a special
# flag to support AVX1.
#
# The following variables are set:
#
# ::
#
#    AVX1_C_FLAGS - flags to add to the C compiler for AVX1 support
#    AVX1_FOUND - true if AVX1 is detected
#
#=============================================================================

set(_AVX1_REQUIRED_VARS)
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ${AVX1_FIND_QUIETLY})

# sample AVX1 source code to test
set(AVX1_C_TEST_SOURCE
"
#include <immintrin.h>
void parasail_memset___m256i(__m256i *b, __m256i c, size_t len)
{
    size_t i;
    for (i=0; i<len; ++i) {
        _mm256_store_si256(&b[i], c);
    }
}
int foo() {
    __m256i vOne = _mm256_set1_epi16(1);
    return _mm256_extract_epi16(vOne, 0);
}
int main(void) { return (int)foo(); }
")

# if these are set then do not try to find them again,
# by avoiding any try_compiles for the flags
if((DEFINED AVX1_C_FLAGS) OR (DEFINED HAVE_AVX1))
else()
  if(WIN32)
    # MSVC can compile AVX intrinsics without the arch flag, however it
    # will detect that AVX code is found and "consider using /arch:AVX".
    set(AVX1_C_FLAG_CANDIDATES
      "/arch:AVX"
    )
  else()
    set(AVX1_C_FLAG_CANDIDATES
      #Empty, if compiler automatically accepts AVX1
      " "
      #clang
      "-mavx"
      #GNU, Intel
      "-march=core-avx-i"
    )
  endif()

  include(CheckCSourceCompiles)

  foreach(FLAG IN LISTS AVX1_C_FLAG_CANDIDATES)
    set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    set(CMAKE_REQUIRED_FLAGS "${FLAG}")
    unset(HAVE_AVX1 CACHE)
    if(NOT CMAKE_REQUIRED_QUIET)
      message(STATUS "Try AVX1 C flag = [${FLAG}]")
    endif()
    check_c_source_compiles("${AVX1_C_TEST_SOURCE}" HAVE_AVX1)
    set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    if(HAVE_AVX1)
      set(AVX1_C_FLAGS_INTERNAL "${FLAG}")
      break()
    endif()
  endforeach()

  unset(AVX1_C_FLAG_CANDIDATES)
  
  set(AVX1_C_FLAGS "${AVX1_C_FLAGS_INTERNAL}"
    CACHE STRING "C compiler flags for AVX1 intrinsics")
endif()

list(APPEND _AVX1_REQUIRED_VARS AVX1_C_FLAGS)

set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})

if(_AVX1_REQUIRED_VARS)
  include(FindPackageHandleStandardArgs)

  find_package_handle_standard_args(AVX1
                                    REQUIRED_VARS ${_AVX1_REQUIRED_VARS})

  mark_as_advanced(${_AVX1_REQUIRED_VARS})

  unset(_AVX1_REQUIRED_VARS)
else()
  message(SEND_ERROR "FindAVX1 requires C or CXX language to be enabled")
endif()
