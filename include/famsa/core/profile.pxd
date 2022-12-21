from libcpp.vector cimport vector

from famsa.core cimport score_t
from famsa.core.params cimport CParams
from famsa.core.sequence cimport CGappedSequence

cdef extern from "core/profile.h" nogil:

    cdef cppclass CProfile:
        CParams* params
        vector[CGappedSequence*] data;
        size_t width
        score_t total_score

        CProfile(CParams* params)
