# coding: utf-8
# cython: language_level=3, linetrace=True
"""Bindings to FAMSA, an algorithm for fast multiple sequence alignments.

References:
    - Deorowicz, S., Debudaj-Grabysz, A., Gudyś, A. (2016)
      *FAMSA: Fast and accurate multiple sequence alignment of huge protein
      families*. Scientific Reports, 6, 33964. :doi:`10.1038/srep33964`.

"""

# --- C imports --------------------------------------------------------------

from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_FORMAT, PyBUF_READ
from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AS_STRING

from libc.stdint cimport uint32_t
from libcpp cimport bool
from libcpp.memory cimport shared_ptr
from libcpp.utility cimport move
from libcpp.vector cimport vector
from libcpp.string cimport string

from famsa.core cimport symbol_t, GAP, GUARD
from famsa.core.params cimport CParams
from famsa.core.sequence cimport CSequence, CGappedSequence
from famsa.msa cimport CFAMSA
from famsa.tree cimport GT, node_t
from famsa.tree.guide_tree cimport GuideTree
from famsa.tree.newick_parser cimport NewickParser
from famsa.tree.abstract_tree_generator cimport AbstractTreeGenerator
from famsa.utils.memory_monotonic cimport memory_monotonic_safe
# from famsa.utils.log cimport Log, LEVEL_NORMAL, LEVEL_DEBUG, LEVEL_VERBOSE


# --- Python imports ---------------------------------------------------------

import os


# --- Constants --------------------------------------------------------------

cdef memory_monotonic_safe* MMA = new memory_monotonic_safe()

cdef char SYMBOLS[25]
for i, x in enumerate(b"ARNDCQEGHILKMFPSTWYVBZX*"):
    SYMBOLS[i] = x

# Log.getInstance(LEVEL_NORMAL).enable()
# Log.getInstance(LEVEL_VERBOSE).enable()
# Log.getInstance(LEVEL_DEBUG).enable()

# --- Utils ------------------------------------------------------------------

cdef extern from *:
    """
    void sort_sequences(vector<CSequence>& sequences) {
        std::stable_sort(sequences.begin(), sequences.end(), [](const CSequence& a, const CSequence& b) -> bool {
          return a.length > b.length || (a.length == b.length && std::lexicographical_compare(a.data, a.data + a.data_size, b.data, b.data + b.data_size));
        });
    }
    void shuffle_sequences(vector<CSequence>& sequences, int shuffle) {
        std::mt19937 mt(shuffle);
        std::shuffle(sequences.begin(), sequences.end(), mt);
    }
    """
    void sort_sequences(vector[CSequence]& sequences)
    void shuffle_sequences(vector[CSequence]& sequences, int shuffle)


# --- Classes ----------------------------------------------------------------

cdef class Sequence:
    """A digitized sequence.
    """

    # --- Magic methods ------------------------------------------------------

    def __init__(self, bytes id, bytes sequence):
        self._cseq = move(CSequence(id, sequence, MMA))
        self._shape[0] = self._cseq.length
        assert self._cseq.mma is not NULL

    def __copy__(self):
        return self.copy()

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        if flags & PyBUF_FORMAT:
            buffer.format = b"b"
        else:
            buffer.format = NULL
        buffer.buf = self._cseq.data
        buffer.internal = NULL
        buffer.itemsize = sizeof(symbol_t)
        buffer.len = self._shape[0] * sizeof(symbol_t)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 1
        buffer.shape = self._shape
        buffer.suboffsets = NULL
        buffer.strides = NULL

    # --- Properties ---------------------------------------------------------

    @property
    def id(self):
        return <bytes> self._cseq.id

    @property
    def sequence(self):
        # code from `CSequence::DecodeSequence`
        cdef uint32_t i
        cdef bytes    seq = PyBytes_FromStringAndSize(NULL, self._cseq.length)
        cdef char*    mem = PyBytes_AS_STRING(seq)

        with nogil:
            for i in range(self._cseq.length):
                if self._cseq.data[i] == GUARD:
                    continue
                elif self._cseq.data[i] == GAP:
                    mem[0] = b'-'
                else:
                    mem[0] = SYMBOLS[self._cseq.data[i]]
                mem = &mem[1]

        return seq

    # --- Methods ------------------------------------------------------------

    cpdef Sequence copy(self):
        """copy(self)\n--

        Copy the sequence data, and return the copy.

        """
        cdef Sequence seq = Sequence.__new__(Sequence)
        seq._cseq = move(CSequence(self._cseq))
        seq._shape = self._shape
        return seq


cdef class GappedSequence:
    """A gapped sequence, storing a single row in an alignment.
    """

    # --- Properties ---------------------------------------------------------

    @property
    def id(self):
        return <bytes> self._gseq.id

    @property
    def sequence(self):
        return <bytes> self._gseq.Decode()


cdef class Alignment:
    """An alignment, stored as a list of `GappedSequence` objects.
    """

    # --- Magic methods ------------------------------------------------------

    def __init__(self):
        self._msa.clear() # create an empty alignment

    def __len__(self):
        return self._msa.size()

    def __getitem__(self, ssize_t index):
        cdef GappedSequence gapped
        cdef ssize_t        index_ = index

        if index_ < 0:
            index_ += self._msa.size()
        if index_ < 0 or index_ >= <ssize_t> self._msa.size():
            raise IndexError(index)

        gapped = GappedSequence.__new__(GappedSequence)
        gapped.alignment = self
        gapped._gseq = self._msa[index_]
        return gapped


cdef class Aligner:

    # --- Magic methods ------------------------------------------------------

    def __cinit__(self):
        self._params = CParams()
        self._params.verbose_mode = True
        self._params.very_verbose_mode = True
        self._params.n_threads = 1

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
            self._params.n_threads = os.cpu_count() or 1
        elif threads >= 1:
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

    # --- Methods ------------------------------------------------------------

    cpdef Alignment align(self, object sequences):
        """align(self, sequences)\n--

        Align sequences together.

        Arguments:
            sequences (iterable of `~pyfamsa.Sequence`): An iterable
                yielding the digitized sequences to align.

        Returns:
            `~pyfamsa.Alignment`: The aligned sequences, in aligned format.

        """
        cdef Sequence                 sequence
        cdef vector[CSequence]        seqvec
        cdef vector[CGappedSequence*] gapvec
        cdef Alignment                alignment = Alignment.__new__(Alignment)

        # create a new aligner
        alignment._famsa = shared_ptr[CFAMSA](new CFAMSA(self._params))
        # copy the aligner input
        for sequence in sequences:
            seqvec.push_back(CSequence(sequence._cseq))

        # align the input and extract the resulting alignment
        with nogil:
            alignment._famsa.get().ComputeMSA(seqvec)
            alignment._famsa.get().GetAlignment(alignment._msa)

        return alignment

    cpdef Tree build_tree(self, object sequences):
        cdef size_t                            i
        cdef Sequence                          sequence
        cdef vector[CSequence]                 seqvec
        cdef vector[CSequence]                 namevec
        cdef vector[CGappedSequence*]          gapvec
        cdef shared_ptr[AbstractTreeGenerator] gen
        cdef CFAMSA*                           famsa    = new CFAMSA(self._params)
        cdef Tree                              tree     = Tree.__new__(Tree)

        # copy the aligner input
        for sequence in sequences:
            seqvec.push_back(CSequence(sequence._cseq))
            tree._names.push_back(CSequence(sequence._cseq.id, string(), NULL))

        # sort or shuffle sequences
        if self._params.shuffle == -1:
            sort_sequences(seqvec)
        else:
            shuffle_sequences(seqvec, self._params.shuffle)

        # set sequence identifiers
        for i in range(seqvec.size()):
            seqvec[i].sequence_no = i

        # generate tree
        try:
            with nogil:
                gen = famsa.createTreeGenerator(self._params)
                gen.get().call( seqvec, tree._tree.raw() )
        finally:
            del famsa

        # return the result tree
        return tree


cdef class Tree:

    def __len__(self):
        return self._tree.raw().size()

    def __getitem__(self, size_t index):
        cdef GappedSequence gapped
        cdef node_t         node
        cdef ssize_t        index_ = index

        if index_ < 0:
            index_ += self._tree.raw().size()
        if index_ < 0 or index_ >= <ssize_t> self._tree.raw().size():
            raise IndexError(index)

        node = self._tree.raw()[index]
        return node

    cpdef bytes dumps(self):
        """dumps(self)\n--

        Dump the tree in Newick format.

        """
        cdef string       out
        cdef NewickParser nw

        with nogil:
            nw.store(self._names, self._tree.raw(), out)

        return <bytes> out
