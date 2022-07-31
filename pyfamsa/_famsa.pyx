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

from libcpp cimport bool
from libcpp.memory cimport shared_ptr
from libcpp.utility cimport move
from libcpp.vector cimport vector
from libcpp.string cimport string

from famsa.core cimport symbol_t
from famsa.core.params cimport CParams
from famsa.core.sequence cimport CSequence, CGappedSequence
from famsa.msa cimport CFAMSA
from famsa.tree cimport GT
from famsa.utils.memory_monotonic cimport memory_monotonic_safe
# from famsa.utils.log cimport Log, LEVEL_NORMAL, LEVEL_DEBUG, LEVEL_VERBOSE


# --- Python imports ---------------------------------------------------------

import os


# --- Constants --------------------------------------------------------------

cdef memory_monotonic_safe* MMA = new memory_monotonic_safe()

# Log.getInstance(LEVEL_NORMAL).enable()
# Log.getInstance(LEVEL_VERBOSE).enable()
# Log.getInstance(LEVEL_DEBUG).enable()

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
        return <bytes> self._cseq.DecodeSequence()

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
        if index_ < 0 or index_ >= self._msa.size():
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
