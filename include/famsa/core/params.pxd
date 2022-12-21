from libc.stdint cimport int64_t, uint32_t, uint64_t
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from famsa.core cimport score_t, instruction_set_t
from famsa.tree cimport GT

cdef extern from "core/params.h" namespace "Refinement" nogil:

    cdef cppclass Mode:
        pass


cdef extern from "core/params.h" namespace "Refinement::Mode" nogil:

    const Mode ON
    const Mode OFF
    const Mode AUTO


cdef extern from "core/params.h" nogil:

    cdef cppclass Refinement:
        string toString(Mo)

    cdef enum Distance:
        indel_div_lcs
        sqrt_indel_div_lcs
        neg_lcs_div_indel
        neg_lcs_div_minlen
        neg_lcs_div_len_corrected

    cdef cppclass CParams:
        score_t gap_open
        score_t gap_ext
        score_t gap_term_open
        score_t gap_term_ext

        uint32_t scaler_div
        uint32_t scaler_log
        int guided_alignment_radius

        bool enable_gap_rescaling
        bool enable_gap_optimization
        bool enable_total_score_calculation

        Mode refinement_mode
        uint32_t n_refinements
        uint32_t thr_refinement
        uint32_t thr_internal_refinement

        GT.Method gt_method
        GT.Heuristic gt_heuristic
        Distance distance
        int heuristic_threshold

        int guide_tree_seed
        int subtree_size
        int sample_size
        float cluster_fraction
        int cluster_iters

        string guide_tree_in_file
        bool export_distances
        bool export_tree
        bool generate_square_matrix
        bool calculate_pid
        bool keepDuplicates

        bool test_ref_sequences
        uint64_t ref_seq_subtree_size
        string ref_file_name

        int64_t shuffle
        uint32_t n_threads

        bool gzippd_output
        int gzip_level

        instruction_set_t instruction_set

        bool verbose_mode
        bool very_verbose_mode

        string input_file_name
        string output_file_name

        vector[vector[score_t]] score_matrix
        vector[score_t] score_vector

        CParams() except +
        bool parse(int argc, char** argv, bool& showExpert) except +
        void show_usage(bool expert) except +
