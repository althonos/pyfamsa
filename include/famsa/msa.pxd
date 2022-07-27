from libcpp cimport bool
from libcpp.vector cimport vector

from famsa.core.param cimport CParams
from famsa.core.sequence cimport
from famsa.tree cimport tree_structure

cdef extern from "msa.h" nogil:

    cdef cppclass CFAMSA:
        CFAMSA(CParams& params) except +

        bool ComputeAlignment(tree_structure& guide_tree) except +
        # bool RefineAlignment(CProfile *&profile_to_refine)
        bool GetAlignment(vector<CGappedSequence*>& result) except +
        score_t GetScore()

        const Statistics& getStatistics()
        Statistics& getStatistics()

        bool ComputeMSA(vector[CSequence]& sequences)
