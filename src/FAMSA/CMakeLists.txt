set(FAMSA_SOURCES

    libs/refresh/active_thread_pool/lib/active_thread_pool.h
    libs/refresh/active_thread_pool/lib/utils.h

    src/msa.cpp
    src/msa_refinement.cpp

    src/tree/AbstractTreeGenerator.cpp
    src/tree/Clustering.cpp
    src/tree/DistanceCalculator.cpp
    src/tree/FastTree.cpp
    src/tree/GuideTree.cpp
    src/tree/MSTPrim.cpp
    src/tree/NeighborJoining.cpp
    src/tree/NewickParser.cpp
    src/tree/SingleLinkage.cpp
    src/tree/UPGMA.cpp

    src/utils/timer.cpp
    src/utils/log.cpp
    src/utils/utils.cpp

    src/core/params.cpp
    src/core/profile.cpp
    # src/core/profile_par.cpp
    src/core/profile_seq.cpp
    src/core/sequence.cpp
    src/core/queues.cpp

    src/core/defs.h
    src/core/params.h
    src/core/profile.h
    src/core/queues.h
    src/core/sequence.h
    src/core/version.h
    src/core/scoring_matrix.h

    src/lcs/lcsbp_classic.h
    src/lcs/lcsbp.h

    src/tree/AbstractTreeGenerator.h
    src/tree/AbstractTreeGenerator.hpp
    src/tree/Chained.h
    src/tree/Clustering.h
    src/tree/DistanceCalculator.h
    src/tree/FastTree.h
    src/tree/GuideTree.h
    src/tree/IPartialGenerator.h
    src/tree/MSTPrim.h
    src/tree/NeighborJoining.h
    src/tree/NewickParser.h
    src/tree/SingleLinkage.h
    src/tree/SingleLinkageQueue.h
    src/tree/TreeDefs.h
    src/tree/UPGMA.h

    src/utils/array.h
    src/utils/conversion.h
    src/utils/cpuid.h
    src/utils/deterministic_random.h
    src/utils/log.h
    src/utils/memory_monotonic.h
    src/utils/meta_oper.h
    src/utils/statistics.h
    src/utils/timer.h
    src/utils/utils.h

    src/msa.h
)

foreach(_file IN ITEMS ${FAMSA_SOURCES})
    cmake_path(GET _file FILENAME _name)
    if(EXISTS ${PROJECT_SOURCE_DIR}/patches/${_name}.patch)
        add_custom_command(
            OUTPUT
                ${_file}
            COMMENT
                "Patching ${_file}"
            COMMAND
                ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/src/scripts/apply_patch.py
                    --input ${PROJECT_SOURCE_DIR}/vendor/FAMSA/${_file}
                    --patch ${PROJECT_SOURCE_DIR}/patches/${_name}.patch
                    --output ${CMAKE_CURRENT_BINARY_DIR}/${_file}
            DEPENDS
                ${PROJECT_SOURCE_DIR}/vendor/FAMSA/${_file}
                ${PROJECT_SOURCE_DIR}/patches/${_name}.patch
        )
    else()
        add_custom_command(
            OUTPUT
                ${_file}
            COMMENT
                "Copying ${_file}"
            COMMAND
                cmake -E copy ${PROJECT_SOURCE_DIR}/vendor/FAMSA/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file}
            DEPENDS
                ${PROJECT_SOURCE_DIR}/vendor/FAMSA/${_file}
        )
    endif()
    set(FAMSA_PATCHED_SOURCES ${FAMSA_PATCHED_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
endforeach()

add_library(famsa STATIC ${FAMSA_PATCHED_SOURCES})
target_include_directories(famsa PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/src)
target_include_directories(famsa PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/libs)
target_include_directories(famsa PUBLIC ${PROJECT_SOURCE_DIR}/vendor/FAMSA/libs) # FIXME?

if(HAVE_AVX2)
    message(STATUS "Building FAMSA with AVX2 support")
    add_library(
        famsa_lcs
        STATIC
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp_classic.cpp
    )
    target_compile_definitions(famsa_lcs PRIVATE -DSIMD=2)
    string(REPLACE " " ";" IMPL_FLAGS ${AVX2_C_FLAGS})
    foreach(_flag IN LISTS IMPL_FLAGS)
        target_compile_options(famsa_lcs PRIVATE ${_flag})
    endforeach()
elseif(HAVE_AVX1)
    message(STATUS "Building FAMSA with AVX1 support")
    add_library(
        famsa_lcs
        STATIC
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp_classic.cpp
    )
    target_compile_definitions(famsa_lcs PRIVATE -DSIMD=1)
    string(REPLACE " " ";" IMPL_FLAGS ${AVX1_C_FLAGS})
    foreach(_flag IN LISTS IMPL_FLAGS)
        target_compile_options(famsa_lcs PRIVATE ${_flag})
    endforeach()
elseif(HAVE_SSE4)
    message(STATUS "Building FAMSA with SSE4 support")
    add_library(
        famsa_lcs
        STATIC
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp_classic.cpp
    )
    target_compile_definitions(famsa_lcs PRIVATE -DSIMD=0)
    string(REPLACE " " ";" IMPL_FLAGS ${SSE4_C_FLAGS})
    foreach(_flag IN LISTS IMPL_FLAGS)
        target_compile_options(famsa_lcs PRIVATE ${_flag})
    endforeach()
elseif(HAVE_NEON)
    message(STATUS "Building FAMSA with NEON support")
    add_library(
        famsa_lcs
        STATIC
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp_classic.cpp
    )
    target_compile_definitions(famsa_lcs PRIVATE -DSIMD=4)
    string(REPLACE " " ";" IMPL_FLAGS ${NEON_C_FLAGS})
    foreach(_flag IN LISTS IMPL_FLAGS)
        target_compile_options(famsa_lcs PRIVATE ${_flag})
    endforeach()
else()
    message(STATUS "Building FAMSA without SIMD code")
    add_library(
        famsa_lcs
        STATIC
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp_classic.cpp
    )
    target_compile_definitions(famsa_lcs PRIVATE -DSIMD=0)
endif()

if(HAVE_POPCNT)
    string(REPLACE " " ";" IMPL_FLAGS ${POPCNT_C_FLAGS})
    foreach(_flag IN LISTS IMPL_FLAGS)
        target_compile_options(famsa_lcs PRIVATE ${_flag})
    endforeach()
endif()

target_link_libraries(famsa PUBLIC famsa_lcs)
target_compile_definitions(famsa PUBLIC _HAS_STD_BYTE=0)  # needed for Windows
