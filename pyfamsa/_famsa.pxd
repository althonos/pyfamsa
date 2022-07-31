# distutils: language = c++
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libcpp.vector cimport vector

from famsa.msa cimport CFAMSA
from famsa.core.params cimport CParams
from famsa.core.sequence cimport CSequence, CGappedSequence
from famsa.utils.memory_monotonic cimport memory_monotonic_safe

# --- Allocator --------------------------------------------------------------

cdef memory_monotonic_safe* MMA = new memory_monotonic_safe()

# --- Classes ----------------------------------------------------------------

cdef class Sequence:
    cdef          CSequence  _cseq
    cdef readonly Py_ssize_t _shape[1]

    cpdef Sequence copy(self)


cdef class GappedSequence:
    cdef Alignment        alignment
    cdef CGappedSequence* _gseq


cdef class Alignment:
    cdef CFAMSA*                  _famsa
    cdef vector[CGappedSequence*] _msa


cdef class Aligner:
    cdef CParams _params

    cpdef Alignment align(self, object sequences)
