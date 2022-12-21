from libc.stdint cimport uint32_t
from libcpp cimport bool
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector

from famsa.core cimport bit_vec_t, symbol_t
from famsa.utils.memory_monotonic cimport memory_monotonic_safe


cdef extern from "core/sequence.h" nogil:

    cdef cppclass CSequence:
        uint32_t length
        uint32_t data_size
        symbol_t* data
        bit_vec_t* p_bit_masks
        uint32_t p_bv_len

        int sequence_no
        int original_no
        string id

        memory_monotonic_safe* mma

        vector[bool] uppercase

        CSequence()
        CSequence(const string& id, const string& seq, int sequence_no, memory_monotonic_safe* mma) except +
        CSequence(CSequence&& x)

        void DataResize(uint32_t new_size, symbol_t new_symbol) except +

        void ComputeBitMasks() except +
        void ReleaseBitMasks() except +
        string DecodeSequence() except +

    cdef struct CSequenceView:
        uint32_t length
        uint32_t padding1
        symbol_t* data

    cdef cppclass CGappedSequence:
        symbol_t* symbols
        size_t size
        size_t symbols_size
        size_t gapped_size
        size_t dps_size
        size_t dps_size_div2

        vector[uint32_t] n_gaps
        vector[uint32_t] dps

        string id
        vector[bool] uppercase

        CGappedSequence(CSequence&& sequence) except +
        CGappedSequence(const CGappedSequence&& sequence) except +
        CGappedSequence(CGappedSequence&& _gapped_sequence)

        void InsertGap(uint32_t pos) except +
        void InsertGaps(uint32_t pos, uint32_t n) except +
        void InsertGapsVector(const vector[pair[uint32_t, uint32_t]]& v_gaps) except +

        void RemoveGap(size_t pos) except +
        void RemoveGaps(size_t pos, uint32_t n) except +
        symbol_t GetSymbol(size_t pos)

        void DecodeRaw(symbol_t* seq) except +
        string Decode() except +
        uint32_t NoSymbols() except +

        void InsertFront(symbol_t new_symbol) except +

        void Clear() except +
        void ClearDPS() except +
