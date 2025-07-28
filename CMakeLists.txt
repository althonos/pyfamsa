cmake_minimum_required(VERSION 3.20)
project(${SKBUILD_PROJECT_NAME} VERSION ${SKBUILD_PROJECT_VERSION} LANGUAGES C CXX)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

message(STATUS "Building ${SKBUILD_PROJECT_NAME} for ${CMAKE_SYSTEM_PROCESSOR}")

# --- Setup include folders ----------------------------------------------------

set_property(GLOBAL PROPERTY PYTHON_EXTENSIONS_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(CYTHON_HEADERS_DIR ${CMAKE_CURRENT_LIST_DIR}/include)

# --- Add global C++ flags -----------------------------------------------------

include(CheckCXXCompilerFlag)

check_cxx_compiler_flag(-funroll-loops HAVE_UNROLL_LOOPS)
if(HAVE_UNROLL_LOOPS)
    add_compile_options(-funroll-loops)
endif()

add_definitions(-DOLD_ATOMIC_FLAG -DNO_PROFILE_PAR)

if(WIN32)
    add_definitions(-D_WIN32)
endif()

# --- Add scripts for Cython ---------------------------------------------------

# include("src/scripts/cmake/FindAVX1.cmake")
# include("src/scripts/cmake/FindAVX2.cmake")
# include("src/scripts/cmake/FindAVX512.cmake")
include("src/scripts/cmake/FindNEON.cmake")
include("src/scripts/cmake/FindPOPCNT.cmake")
include("src/scripts/cmake/FindSSE2.cmake")
include("src/scripts/cmake/FindSSE4.cmake")
include("src/scripts/cmake/CythonExtension.cmake")

# --- Compile code -------------------------------------------------------------

add_subdirectory("src")
