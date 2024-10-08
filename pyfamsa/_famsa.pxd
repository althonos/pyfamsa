# distutils: language = c++
# cython: language_level=3, linetrace=True

# --- C imports --------------------------------------------------------------

from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector

from famsa.msa cimport CFAMSA
from famsa.core cimport NO_AMINOACIDS
from famsa.core.params cimport CParams
from famsa.core.sequence cimport CSequence, CGappedSequence
from famsa.tree.guide_tree cimport GuideTree as CGuideTree
from famsa.utils.memory_monotonic cimport memory_monotonic_safe

from scoring_matrices.lib cimport ScoringMatrix

# --- Allocator --------------------------------------------------------------

cdef memory_monotonic_safe* MMA = new memory_monotonic_safe()

# --- Classes ----------------------------------------------------------------

cdef class Sequence:
    cdef          shared_ptr[CSequence]  _cseq
    cdef readonly Py_ssize_t             _shape[1]

    cpdef Sequence copy(self)


cdef class GappedSequence:
    cdef shared_ptr[CGappedSequence] _gseq

    cpdef GappedSequence copy(self)


cdef class Alignment:
    cdef vector[shared_ptr[CGappedSequence]] _msa

    cpdef Alignment copy(self)


cdef class Aligner:
    cdef          CParams       _params
    cdef readonly ScoringMatrix scoring_matrix

    cdef int _copy_matrix(self, CFAMSA* famsa) except 1 nogil

    cpdef Alignment align(self, object sequences)
    cpdef Alignment align_profiles(self, Alignment profile1, Alignment profile2)
    cpdef GuideTree build_tree(self, object sequences)


cdef class GuideTree:
    cdef CGuideTree        _tree
    cdef vector[CSequence] _names

    cpdef bytes dumps(self)
    cpdef ssize_t dump(self, object file) except -1
