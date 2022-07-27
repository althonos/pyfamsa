from libc.stdint cimport int64_t

cdef extern from "core/defs.h" nogil:

    ctypedef int64_t score_t

    cdef enum instruction_set_t:
        none
        sse
        sse2
        sse3
        sse3s
        sse41
        sse42
        avx
        avx2

    ctypedef char symbol_t
    ctypedef int counter_t

    ctypedef unsigned long long bit_vec_t

    const int bv_size
    const int bv_size128
    const int bv_size256

    const score_t infty
    const double cost_cast_factor

    const symbol_t GAP
    const symbol_t GAP_OPEN
    const symbol_t GAP_EXT
    const symbol_t GAP_TERM_EXT
    const symbol_t GAP_TERM_OPEN
    const symbol_t UNKNOWN_SYMBOL

    const size_t NO_SYMBOLS
    const symbol_t GUARD
    const symbol_t NO_AMINOACIDS
    const symbol_t NO_VALID_AMINOACIDS
    const symbol_t NO_AMINOACIDS_AND_GAPS
    const symbol_t NO_AA_SYMBOLS
