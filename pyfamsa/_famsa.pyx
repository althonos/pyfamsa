# coding: utf-8
# cython: language_level=3, linetrace=True

from libcpp cimport bool
from libcpp.utility cimport move
from libcpp.vector cimport vector
from libcpp.string cimport string

from famsa.msa cimport CFAMSA
from famsa.core.io_service cimport IOService
from famsa.core.params cimport CParams
from famsa.core.sequence cimport CSequence, CGappedSequence
from famsa.tree cimport GT
from famsa.utils.memory_monotonic cimport memory_monotonic_safe
from famsa.utils.log cimport Log, LEVEL_NORMAL, LEVEL_DEBUG, LEVEL_VERBOSE

import os

cdef memory_monotonic_safe* MMA = new memory_monotonic_safe()


cdef class Sequence:

    cdef CSequence _cseq

    def __init__(self, bytes id, bytes sequence):
        self._cseq = move(CSequence(id, sequence, MMA))
        assert self._cseq.mma is not NULL

    def __copy__(self):
        return self.copy()

    @property
    def id(self):
        return <bytes> self._cseq.id

    @property
    def sequence(self):
        return <bytes> self._gseq.DecodeSequence()

    cpdef Sequence copy(self):
        cdef Sequence seq = Sequence.__new__(Sequence)
        seq._cseq = move(CSequence(self._cseq))
        return seq


cdef class GappedSequence:

    cdef Alignment        alignment
    cdef CGappedSequence* _gseq

    @property
    def id(self):
        return <bytes> self._gseq.id

    @property
    def sequence(self):
        return <bytes> self._gseq.Decode()


cdef class Alignment:

    cdef CFAMSA*                  _famsa
    cdef vector[CGappedSequence*] _msa

    def __len__(self):
        return self._msa.size()

    def __getitem__(self, ssize_t index):
        cdef GappedSequence gapped
        cdef ssize_t        index_ = index

        if index_ < 0:
            index_ += self._msa.size()
        if index_ < 0 or index_ >= self._msa.size():
            raise IndexError(index)

        gapped = GappedSequence.__new__(GappedSequence)
        gapped.alignment = self
        gapped._gseq = self._msa[index_]
        return gapped


cdef class Aligner:

    cdef CParams                _params

    def __cinit__(self):
        self._params = CParams()
        self._params.verbose_mode = True
        self._params.very_verbose_mode = True
        self._params.n_threads = 1

    def __dealloc__(self):
        # del self._mma
        pass

    def __init__(
        self,
        *,
        int threads=0,
        object guide_tree="sl",
        object tree_heuristic=None,
        int medoid_threshold=0,
        int n_refinements=100,
        bool force_refinement=False,
    ):
        """__init__(self, *, threads=0, guide_tree="sl", tree_heuristic=None, medoid_threshold=0, n_refinements=100, force_refinement=False)\n--

        Create a new aligner with the given configuration.

        Keyword Arguments:
            threads (`int`): The number of threads to use for parallel
                computations. If *0* given (the default), use `os.cpu_count`
                to spawn one thread per CPU on the host machine.
            guide_tree (`str`): The method for building the guide tree.
                Supported values are: ``sl`` for MST+Prim single linkage,
                ``slink`` for SLINK single linkage, ``upgma`` for UPGMA,
                ``nj`` for neighbour joining.
            tree_heuristic (`str` or `None`): The heuristic to use for
                constructing the tree. Supported values are: ``medoid`` for
                medoid trees, ``part`` for part trees, or `None` to disable
                heuristics.
            medoid_threshold (`int`): The minimum number of sequences a
                set must contain for medoid trees to be used, if enabled
                with ``tree_heuristic``.
            n_refinements (`int`): The number of refinement iterations to
                run.
            force_refinement (`bool`): Set to `True` to force refinement
                on sequence sets larger than 1000 sequences.

        """
        if threads == 0:
            self._params.n_threads = os.cpu_count()
        elif threads > 1:
            self._params.n_threads = threads
        else:
            raise ValueError("`threads` argument must be positive")

        if guide_tree == "sl":
            self._params.gt_method = GT.Method.MST_Prim
        elif guide_tree == "slink":
            self._params.gt_method = GT.Method.SLINK
        elif guide_tree == "upgma":
            self._params.gt_method = GT.Method.UPGMA
        elif guide_tree == "nj":
            self._params.gt_method = GT.Method.NJ
        else:
            raise ValueError(f"Invalid value for `guide_tree` argument: {guide_tree!r}")

        if tree_heuristic is None:
            self._params.gt_heuristic = GT.Heuristic.None
        elif tree_heuristic == 'medoid':
            self._params.gt_heuristic = GT.Heuristic.ClusterTree
        elif tree_heuristic == 'parttree':
            self._params.gt_heuristic = GT.Heuristic.PartTree
        else:
            raise ValueError(f"Invalid value for `tree_heuristic` argument: {tree_heuristic!r})")

        if medoid_threshold >= 0:
            self._params.heuristic_threshold = medoid_threshold
        else:
            raise ValueError("`medoid_threshold` argument must be positive")

        if n_refinements >= 0:
            self._params.n_refinements = n_refinements
        else:
            raise ValueError("`n_refinements` argument must be positive")

        if force_refinement:
            self._params.enable_auto_refinement = False

    cpdef Alignment align(self, object sequences):

        cdef Sequence                 sequence
        cdef vector[CSequence]        seqvec
        cdef vector[CGappedSequence*] gapvec
        cdef Alignment                alignment = Alignment.__new__(Alignment)

        alignment._famsa = new CFAMSA(self._params)

        for sequence in sequences:
            seqvec.push_back(CSequence(sequence._cseq))

        alignment._famsa.ComputeMSA(seqvec)
        alignment._famsa.GetAlignment(alignment._msa)

        return alignment

# Log.getInstance(LEVEL_NORMAL).enable()
# Log.getInstance(LEVEL_VERBOSE).enable()
# Log.getInstance(LEVEL_DEBUG).enable()


# def align(object input_file, object output_file):
#     cdef string in_  = os.fsencode(input_file)
#     cdef string out_ = os.fsencode(output_file)
#
#     cdef vector[CSequence] sequences
#     print(IOService.loadFasta(in_, sequences))
#
#
#     cdef CParams params   = CParams()
#     params.n_threads = 8
#
#     cdef CFAMSA* famsa    = new CFAMSA(params)
#     famsa.ComputeMSA(sequences)
#
#     cdef vector[CGappedSequence*] gapvec
#     famsa.GetAlignment(gapvec)
#
#     IOService.saveAlignment(out_, gapvec, True, -1)
