from libcpp.vector cimport vector

from famsa.core cimport score_t
from famsa.core.params cimport CParams
from famsa.core.sequence cimport CGappedSequence

cdef extern from "core/scoring_matrix.h" namespace "ScoringMatrices" nogil:

    ctypedef int matrix_type_t
        # pass

cdef extern from "core/scoring_matrix.h" namespace "ScoringMatrices::matrix_type_t" nogil:

    const matrix_type_t MIQS
    const matrix_type_t PFASUM31
    const matrix_type_t PFASUM43
    const matrix_type_t PFASUM60
        

cdef extern from "core/scoring_matrix.h" nogil:

    cppclass ScoringMatrices:
        @staticmethod
        double[24][24]& get_matrix(const matrix_type_t matrix_type)

    # cdef cppclass CProfile:
    #     CParams* params
    #     vector[CGappedSequence*] data
    #     size_t width
    #     score_t total_score

    #     CProfile(CParams* params)
