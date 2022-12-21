from libcpp cimport bool
from libcpp.memory cimport shared_ptr
from libcpp.vector cimport vector

from famsa.core cimport score_t, instruction_set_t
from famsa.core.params cimport CParams
from famsa.core.profile cimport CProfile
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

        void RefineRandom(CProfile* profile_to_refine, vector[size_t]& dest_prof_id) except +
        void RefineMostEmptyAndFullColumn(CProfile* profile_to_refine, vector[size_t]& dest_prof_id, vector[size_t]& gap_stats, bool valid_gap_stats) except +

        shared_ptr[AbstractTreeGenerator] createTreeGenerator(const CParams& params) except +

        void sortAndExtendSequences(vector[CSequence]& sequences) except +
        void extendSequences(vector[CSequence]& sequences) except +
        void shrinkSequences(vector[CSequence]& sequences) except +
        void removeDuplicates(vector[CSequence*]& sorted_seqs, vector[int]& original2sorted) except +

        bool ComputeMSA(vector[CSequence]& sequences) except +

