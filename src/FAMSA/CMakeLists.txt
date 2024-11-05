set(FAMSA_SOURCES
    msa.cpp
    msa_refinement.cpp

    tree/AbstractTreeGenerator.cpp
    tree/Clustering.cpp
    tree/DistanceCalculator.cpp
    tree/FastTree.cpp
    tree/GuideTree.cpp
    tree/MSTPrim.cpp
    tree/NeighborJoining.cpp
    tree/NewickParser.cpp
    tree/SingleLinkage.cpp
    tree/UPGMA.cpp

    utils/timer.cpp
    utils/log.cpp
    utils/utils.cpp

    core/params.cpp
    core/profile.cpp
    core/profile_par.cpp
    core/profile_seq.cpp
    core/sequence.cpp
    core/queues.cpp

    core/defs.h
    core/params.h
    core/profile.h
    core/queues.h
    core/sequence.h
    core/version.h

    lcs/lcsbp_avx_intr.h
    lcs/lcsbp_avx2_intr.h
    lcs/lcsbp_classic.h
    lcs/lcsbp_neon_intr.h
    lcs/lcsbp.h

    tree/AbstractTreeGenerator.h
    tree/AbstractTreeGenerator.hpp
    tree/Chained.h
    tree/Clustering.h
    tree/DistanceCalculator.h
    tree/FastTree.h
    tree/GuideTree.h
    tree/IPartialGenerator.h
    tree/MSTPrim.h
    tree/NeighborJoining.h
    tree/NewickParser.h
    tree/SingleLinkage.h
    tree/SingleLinkageQueue.h
    tree/TreeDefs.h
    tree/UPGMA.h

    utils/array.h
    utils/conversion.h
    utils/cpuid.h
    utils/deterministic_random.h
    utils/log.h
    utils/memory_monotonic.h
    utils/meta_oper.h
    utils/pooled_threads.h
    utils/statistics.h
    utils/timer.h
    utils/utils.h

    msa.h
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
                    --input ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/${_file}
                    --patch ${PROJECT_SOURCE_DIR}/patches/${_name}.patch
                    --output ${CMAKE_CURRENT_BINARY_DIR}/${_file}
            DEPENDS
                ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/${_file}
                ${PROJECT_SOURCE_DIR}/patches/${_name}.patch
        )
    else()
        add_custom_command(
            OUTPUT
                ${_file}
            COMMAND
                cmake -E copy ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file}
            DEPENDS
                ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/${_file}
        )
    endif()
    set(FAMSA_PATCHED_SOURCES ${FAMSA_PATCHED_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
endforeach()

add_library(famsa STATIC ${FAMSA_PATCHED_SOURCES})
target_include_directories(famsa PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(famsa PUBLIC ${PROJECT_SOURCE_DIR}/vendor/FAMSA/libs) # FIXME?

if(HAVE_AVX2)
    message(STATUS "Building FAMSA with AVX2 support")
    add_library(
        famsa_lcs
        STATIC
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp_classic.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp_avx_intr.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp_avx2_intr.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/utils/utils_avx.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/utils/utils_avx2.cpp
    )
    target_compile_definitions(famsa_lcs PRIVATE -DSIMD=2)
    target_compile_options(famsa_lcs PRIVATE ${AVX2_C_FLAGS})
elseif(HAVE_AVX1)
    message(STATUS "Building FAMSA with AVX1 support")
    add_library(
        famsa_lcs
        STATIC
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp_classic.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp_avx_intr.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/utils/utils_avx.cpp
    )
    target_compile_definitions(famsa_lcs PRIVATE -DSIMD=1)
    target_compile_options(famsa_lcs PRIVATE ${AVX1_C_FLAGS})
elseif(HAVE_SSE4)
    message(STATUS "Building FAMSA with SSE4 support")
    add_library(
        famsa_lcs
        STATIC
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp_classic.cpp
    )
    target_compile_definitions(famsa_lcs PRIVATE -DSIMD=0)
    target_compile_options(famsa_lcs PRIVATE ${SSE4_C_FLAGS})
elseif(HAVE_NEON)
    message(STATUS "Building FAMSA with NEON support")
    add_library(
        famsa_lcs
        STATIC
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp_classic.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/lcs/lcsbp_neon_intr.cpp
        ${PROJECT_SOURCE_DIR}/vendor/FAMSA/src/utils/utils_neon.cpp
    )
    target_compile_definitions(famsa_lcs PRIVATE -DSIMD=4)
    target_compile_options(famsa_lcs PRIVATE ${NEON_C_FLAGS})
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
    target_compile_options(famsa_lcs PRIVATE ${POPCNT_C_FLAGS})    
endif()


target_link_libraries(famsa PUBLIC famsa_lcs)

