from libcpp cimport bool
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector

from famsa.core cimport score_t, instruction_set_t
from famsa.core.params cimport CParams
from famsa.core.sequence cimport CSequence, CGappedSequence
from famsa.tree cimport tree_structure
from famsa.tree.abstract_tree_generator cimport AbstractTreeGenerator

cdef extern from "msa.h" nogil:

    cdef cppclass CFAMSA:
        instruction_set_t instruction_set
        score_t avg_sim

        CFAMSA(CParams& params) except +

        bool ComputeAlignment(tree_structure& guide_tree) except +
        # bool RefineAlignment(CProfile *&profile_to_refine)
        bool GetAlignment(vector[CGappedSequence*]& result) except +
        score_t GetScore()

        # const Statistics& getStatistics()
        # Statistics& getStatistics()

        bool ComputeMSA(vector[CSequence]& sequences)  except +

        shared_ptr[AbstractTreeGenerator] createTreeGenerator(const CParams& params) except +
