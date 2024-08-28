# coding: utf-8
# cython: language_level=3, linetrace=True
"""Bindings to FAMSA, an algorithm for fast multiple sequence alignments.

Attributes:
    FAMSA_ALPHABET (`str`): The alphabet used by default by FAMSA to encode
        sequences with ordinal encoding.
    MIQS (`~scoring_matrices.ScoringMatrix`): The MIQS scoring matrix proposed
        by Yamada & Tomii (2014), used by default in FAMSA for scoring
        alignments.

References:
    - Deorowicz, S., Debudaj-Grabysz, A., Gudyś, A. (2016)
      *FAMSA: Fast and accurate multiple sequence alignment of huge protein
      families*. Scientific Reports, 6, 33964. :doi:`10.1038/srep33964`.
    - Yamada, K., Tomii, K. (2014).
      *Revisiting amino acid substitution matrices for identifying distantly
      related proteins*. Bioinformatics (Oxford, England), 30(3), 317-325.
      :doi:`10.1093/bioinformatics/btt694`. :pmid:`24281694`.

"""

# --- C imports --------------------------------------------------------------

from cpython cimport Py_buffer
from cpython.memoryview cimport PyMemoryView_FromMemory
from cpython.buffer cimport PyBUF_FORMAT, PyBUF_READ
from cpython.bytes cimport PyBytes_FromStringAndSize, PyBytes_AS_STRING

from libc.stdint cimport uint32_t, SIZE_MAX
from libc.string cimport memset
from libc.math cimport roundf
from libcpp cimport bool
from libcpp.algorithm cimport transform
from libcpp.memory cimport shared_ptr, make_shared
from libcpp.utility cimport move
from libcpp.vector cimport vector
from libcpp.string cimport string

cimport famsa.core.version
from famsa.core cimport score_t, symbol_t, GAP, GUARD, NO_AMINOACIDS, cost_cast_factor
from famsa.core.params cimport CParams, ON, OFF, AUTO
from famsa.core.sequence cimport CSequence, CGappedSequence
from famsa.msa cimport CFAMSA, SM_MIQS
from famsa.tree cimport GT, node_t
from famsa.tree.guide_tree cimport GuideTree as CGuideTree
from famsa.tree.newick_parser cimport NewickParser
from famsa.tree.abstract_tree_generator cimport AbstractTreeGenerator
from famsa.utils.memory_monotonic cimport memory_monotonic_safe
# from famsa.utils.log cimport Log, LEVEL_NORMAL, LEVEL_DEBUG, LEVEL_VERBOSE

from scoring_matrices.lib cimport ScoringMatrix

# --- Python imports ---------------------------------------------------------

import datetime
import collections
import os

include "_version.py"

# --- Constants --------------------------------------------------------------

FAMSA_ALPHABET = "ARNDCQEGHILKMFPSTWYVBZX*"

cdef memory_monotonic_safe* MMA = new memory_monotonic_safe()

cdef char SYMBOLS[NO_AMINOACIDS]
for i, x in enumerate(FAMSA_ALPHABET):
    SYMBOLS[i] = ord(x)

cdef ScoringMatrix _make_miqs():
    cdef list row
    cdef list weights = []
    for i in range(NO_AMINOACIDS):
        row = []
        for j in range(NO_AMINOACIDS):
            row.append(round(SM_MIQS[i][j], 4))
        weights.append(row)
    return ScoringMatrix(
        weights,
        alphabet=FAMSA_ALPHABET,
        name="MIQS",
    )

MIQS = _make_miqs()

# Log.getInstance(LEVEL_NORMAL).enable()
# Log.getInstance(LEVEL_VERBOSE).enable()
# Log.getInstance(LEVEL_DEBUG).enable()


# --- Version info about the wrapped FAMSA code ------------------------------

def famsa_info():
    """Get information about the embedded FAMSA version.

    This function returns an anonymous named tuple with attributes
    ``version``, ``version_info``, ``date`` and ``authors`` that can
    be used to get the build information from FAMSA.

    Example:
        >>> info = famsa_info()
        >>> info.version
        '2.2.2'
        >>> info.date
        datetime.date(2022, 10, 9)

    """
    _VersionInfo = collections.namedtuple("_VersionInfo", ["major", "minor", "micro"])
    _Info = collections.namedtuple("_Info", ["version", "version_info", "date", "authors"])
    major, minor, micro = map(int, famsa.core.version.FAMSA_VER.split(b"."))
    version_info = _VersionInfo(major=major, minor=minor, micro=micro)
    date = datetime.datetime.strptime(famsa.core.version.FAMSA_DATE.decode(), '%Y-%m-%d').date()
    version = famsa.core.version.FAMSA_VER.decode()
    authors = list(map(str.strip, famsa.core.version.FAMSA_AUTHORS.decode().split(",")))
    return _Info(version=version, version_info=version_info, date=date, authors=authors)


# --- Utils ------------------------------------------------------------------

cdef extern from *:
    """
    void sort_sequences(vector<CSequence*>& sequences) {
        std::stable_sort(sequences.begin(), sequences.end(), [](const CSequence* a, const CSequence* b) -> bool {
          return a->length > b->length || (a->length == b->length && std::lexicographical_compare(a->data, a->data + a->data_size, b->data, b->data + b->data_size));
        });
    }
    void shuffle_sequences(vector<CSequence*>& sequences, int shuffle) {
        std::mt19937 mt(shuffle);
        std::shuffle(sequences.begin(), sequences.end(), mt);
    }
    """
    void sort_sequences(vector[CSequence*]& sequences)
    void shuffle_sequences(vector[CSequence*]& sequences, int shuffle)

# --- PyPI patch --------------------------------------------------------------

cdef extern from *:
    """
    #ifndef HAS_PYINTERPRETERSTATE_GETID
    int64_t PyInterpreterState_GetID(PyInterpreterState *interp) {
        return 0;
    }
    #endif
    """

# --- Classes ----------------------------------------------------------------

cdef class Sequence:
    """A digitized sequence.
    """

    # --- Magic methods ------------------------------------------------------

    def __init__(self, bytes id, bytes sequence):
        """__init__(self, id, sequence)\n--\n

        Create a new sequence.

        Arguments:
            id (`bytes`): The sequence identifier.
            sequence (`bytes`): The sequence contents.

        Raises:
            `ValueError`: when initializing an empty sequence.

        """
        if len(sequence) == 0:
            raise ValueError("Cannot create an empty sequence")
        self._cseq.reset(new CSequence(id, sequence, 0, NULL))
        self._shape[0] = self._cseq.get().length

    def __copy__(self):
        return self.copy()

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        if flags & PyBUF_FORMAT:
            buffer.format = b"b"
        else:
            buffer.format = NULL
        buffer.buf = self._cseq.get().data
        buffer.internal = NULL
        buffer.itemsize = sizeof(symbol_t)
        buffer.len = self._shape[0] * sizeof(symbol_t)
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = 1
        buffer.shape = self._shape
        buffer.suboffsets = NULL
        buffer.strides = NULL

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.id}, {self.sequence})"

    def __reduce__(self):
        return type(self), (self.id, self.sequence)

    # --- Properties ---------------------------------------------------------

    @property
    def id(self):
        """`bytes`: The identifier of the sequence.
        """
        return <bytes> self._cseq.get().id

    @property
    def sequence(self):
        """`bytes`: The symbols of the sequence as an ASCII string.
        """
        # code from `CSequence::DecodeSequence`
        cdef uint32_t         i
        cdef bytes            seq 
        cdef char*            mem
        cdef const CSequence* cseq

        cseq = self._cseq.get()
        seq = PyBytes_FromStringAndSize(NULL, cseq.length)
        mem = PyBytes_AS_STRING(seq)

        with nogil:
            for i in range(cseq.length):
                if cseq.data[i] == GUARD:
                    continue
                elif cseq.data[i] == GAP:
                    mem[0] = b'-'
                else:
                    mem[0] = SYMBOLS[cseq.data[i]]
                mem += 1

        return seq

    @property
    def size(self):
        """`int`: The number of symbols in the sequence.
        """
        return self._gseq.get().size

    # --- Methods ------------------------------------------------------------

    cpdef Sequence copy(self):
        """Copy the sequence data, and return the copy.
        """
        cdef Sequence seq = Sequence.__new__(Sequence)
        seq._cseq = self._cseq
        seq._shape = self._shape
        return seq


cdef class GappedSequence:
    """A gapped sequence, storing a single row in an alignment.

    Internally, gapped sequences are not stored in a single string in FAMSA,
    but using a gap counter for each position. This allows saving a lot of
    memory for alignments containing lots of gaps, but adds some overhead to
    decode the full sequence.

    """

    # --- Magic methods ------------------------------------------------------

    def __init__(self, bytes id, bytes sequence):
        """__init__(self, id, sequence)\n--\n

        Create a new sequence.

        Arguments:
            id (`bytes`): The sequence identifier.
            sequence (`bytes`): The sequence contents.

        Raises:
            `ValueError`: when initializing an empty sequence.

        """
        if len(sequence) == 0:
            raise ValueError("Cannot create an empty sequence")
        self._gseq.reset(new CGappedSequence(id, sequence, 0, NULL))

    def __repr__(self):
        cdef str ty = self.__class__.__name__
        return f"{ty}({self.id!r}, {self.sequence!r})"

    def __reduce__(self):
        return type(self), (self.id, self.sequence)

    def __copy__(self):
        return self.copy()

    # --- Properties ---------------------------------------------------------

    @property
    def id(self):
        """`bytes`: The identifier of the gapped sequence.
        """
        return <bytes> self._gseq.get().id

    @property
    def sequence(self):
        """`bytes`: The symbols of the gapped sequence as an ASCII string.
        """
        # code from `CSequence::DecodeSequence` without inplace modification
        cdef uint32_t               i
        cdef char                   symbol
        cdef const CGappedSequence* gseq
        cdef bytes                  seq
        cdef char*                  mem

        gseq = self._gseq.get()
        seq = PyBytes_FromStringAndSize(NULL, gseq.gapped_size)
        mem = PyBytes_AS_STRING(seq)

        with nogil:

            # starting gaps
            memset(mem, b'-', gseq.n_gaps[0])
            mem += gseq.n_gaps[0]

            for i in range(1, gseq.size+1):
                symbol = SYMBOLS[gseq.symbols[i]]
                if not gseq.uppercase[i - 1]:
                    symbol += 32

                mem[0] = symbol
                mem += 1

                memset(mem, b'-', gseq.n_gaps[i])
                mem += gseq.n_gaps[i]

        return seq

    @property
    def size(self):
        """`int`: The number of symbols in the sequence, excluding gaps.
        """
        return self._gseq.get().size

    @property
    def gapped_size(self):
        """`int`: The number of symbols in the sequence, including gaps.
        """
        return self._gseq.get().gapped_size

    # --- Methods ------------------------------------------------------------

    cpdef GappedSequence copy(self):
        """Copy the sequence data, and return the copy.
        """
        cdef GappedSequence gseq = GappedSequence.__new__(GappedSequence)
        gseq._gseq = self._gseq
        return gseq


cdef class Alignment:
    """An alignment, stored as a list of `GappedSequence` objects.

    An alignment can be created from an iterable of `GappedSequence` that
    are all the same gapped size::

        >>> s1 = GappedSequence(b"seq1", b"MA-WMRLLPL")
        >>> s2 = GappedSequence(b"seq2", b"MALWTR-RPL")
        >>> alignment = Alignment([s1, s2])

    The individual rows of the alignment can be accessed using the usual
    indexing operation, and the alignment can be sliced to select only
    certain rows::

        >>> alignment[0].id
        b'seq1'
        >>> len(alignment[:1])
        1

    """

    # --- Magic methods ------------------------------------------------------

    def __init__(self, object sequences = ()):
        """__init__(self, sequences=())\n--\n

        Create a new alignment containing the provided sequences.

        Arguments:
            sequences (iterable of `GappedSequence`): The gapped sequences
                to store in the alignment.

        Raises:
            `ValueError`: When the given sequences are not all of the same
                gapped size, thus are not forming a proper alignment.

        """
        cdef GappedSequence   gseq
        cdef CGappedSequence* seq
        cdef int              i
        cdef size_t           gsize
        cdef size_t           width = SIZE_MAX

        self._msa.clear()
        for i, gseq in enumerate(sequences):
            gsize = gseq._gseq.get().gapped_size
            if width == SIZE_MAX:
                width = gsize
            elif gsize != width:
                raise ValueError(f"sequence sizes mismatch: {gsize} != {width}")
            # seq = new CGappedSequence(gseq._gseq[0])
            # seq.original_no = seq.sequence_no = i
            self._msa.push_back(gseq._gseq)

    def __reduce__(self):
        return type(self), (list(self),)

    def __len__(self):
        return self._msa.size()

    def __getitem__(self, object index):
        cdef GappedSequence   gapped
        cdef CGappedSequence* gseq
        cdef ssize_t          index_
        cdef Alignment        ali
        cdef size_t           length = self._msa.size()

        if isinstance(index, slice):
            indices = range(*index.indices(length))
            ali = Alignment.__new__(Alignment)
            for i in indices:
                ali._msa.push_back(self._msa[i])
            return ali
        else:
            index_ = index
            if index_ < 0:
                index_ += self._msa.size()
            if index_ < 0 or index_ >= <ssize_t> self._msa.size():
                raise IndexError(index)
            gapped = GappedSequence.__new__(GappedSequence)
            gapped._gseq = self._msa[index_]
            return gapped

    # --- Methods ------------------------------------------------------------

    cpdef Alignment copy(self):
        """Copy the sequence object, and return the copy.
        """
        cdef shared_ptr[CGappedSequence] gseq
        cdef Alignment                   ali  = Alignment.__new__(Alignment)

        with nogil:
            for gseq in self._msa:
                ali._msa.push_back(gseq)

        return ali

cdef class Aligner:
    """A single FAMSA aligner.

    Attributes:
        scoring_matrix (`~scoring_matrices.ScoringMatrix`): The scoring
            matrix used for scoring alignments.

    .. versionadded:: 0.4.0
       The ``scoring_matrix`` attribute.

    """

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
        bool keep_duplicates=False,
        object refine=None,
        object scoring_matrix=None,
    ):
        """__init__(self, *, threads=0, guide_tree="sl", tree_heuristic=None, medoid_threshold=0, n_refinements=100, keep_duplicates=False, refine=None, scoring_matrix=None)\n--\n
        
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
            keep_duplicates (`bool`): Set to `True` to avoid discarding
                duplicate sequences before building trees or alignments.
            refine (`bool` or `None`): Set to `True` to force refinement,
                `False` to disable refinement, or leave as `None` to disable
                refinement automatically for sets of more than 1000 sequences.
            scoring_matrix (`~scoring_matrices.ScoringMatrix` or `str`): The
                scoring matrix to use for scoring alignments. By default, the
                *MIQS* matrix by Yamada & Tomii (2014) is used like in the
                original FAMSA implementation.

        .. versionadded:: 0.4.0
           The ``scoring_matrix`` argument.

        """
        self._params.keepDuplicates = keep_duplicates

        if refine is True:
            self._params.refinement_mode = ON
        elif refine is False:
            self._params.refinement_mode = OFF
        elif refine is None:
            self._params.refinement_mode = AUTO
        else:
            raise ValueError(f"Invalid value for `refine` argument: {refine!r}")

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

        if scoring_matrix is None:
            self.scoring_matrix = MIQS
        elif isinstance(scoring_matrix, str):
            self.scoring_matrix = ScoringMatrix.from_name(scoring_matrix).shuffle(FAMSA_ALPHABET)
        elif isinstance(scoring_matrix, ScoringMatrix):
            if scoring_matrix.alphabet != FAMSA_ALPHABET:
                raise ValueError(f"invalid scoring matrix alphabet: expected {FAMSA_ALPHABET!r}, got {scoring_matrix.alphabet!r}")
            self.scoring_matrix = scoring_matrix
        else:
            ty = type(scoring_matrix).__name__
            raise TypeError(f"expected str or ScoringMatrix, found {ty}")

    # --- Methods ------------------------------------------------------------

    cdef int _copy_matrix(self, CFAMSA* famsa) except 1 nogil:
        cdef size_t i
        cdef size_t j
        cdef const float** matrix = self.scoring_matrix.matrix_ptr()
        for i in range(<size_t> NO_AMINOACIDS):
            famsa.score_vector[i] = <score_t> roundf(cost_cast_factor * matrix[i][i])
            for j in range(<size_t> NO_AMINOACIDS):
                famsa.score_matrix[i][j] = <score_t> roundf(cost_cast_factor * matrix[i][j])
        return 0

    cpdef Alignment align(self, object sequences):
        """Align sequences together.

        Arguments:
            sequences (iterable of `~pyfamsa.Sequence`): An iterable
                yielding the digitized sequences to align.

        Returns:
            `~pyfamsa.Alignment`: The aligned sequences, in aligned format.

        """
        cdef int                      i
        cdef int                      j
        cdef CGappedSequence*         aligned
        cdef Sequence                 sequence
        cdef CSequence                cseq
        cdef vector[CSequence]        seqvec
        cdef vector[CGappedSequence*] gapvec
        cdef const float**            matrix
        cdef Alignment                alignment = Alignment.__new__(Alignment)
        cdef CFAMSA*                  famsa     = new CFAMSA(self._params)

        try:
            # copy the aligner input and record sequence order
            for i, sequence in enumerate(sequences):
                cseq = CSequence(sequence._cseq.get()[0])
                cseq.sequence_no = cseq.original_no = i
                seqvec.push_back(move(cseq))
            with nogil:
                # copy score matrix weights
                self._copy_matrix(famsa)
                # align the input and extract the resulting alignment
                if seqvec.size() > 0:
                    if not famsa.ComputeMSA(seqvec):
                        raise RuntimeError("failed to align sequences")
                    # take ownership of the final alignment
                    for aligned in famsa.final_profile.data:
                        alignment._msa.emplace_back(aligned)
                    famsa.final_profile.data.clear()
        finally:
            del famsa

        return alignment

    cpdef Alignment align_profiles(self, Alignment profile1, Alignment profile2):
        """Align two profiles together.

        Profile-profile alignment computes a new alignment using sequences
        from the two input alignments while preserving the columns of each
        profile.

        Arguments:
            profile1 (`~pyfamsa.Alignment`): The first profile to align.
            profile2 (`~pyfamsa.Alignment`): The second profile to align.

        Returns:
            `~pyfamsa.Alignment`: The resulting profile-profile alignment.

        .. versionadded:: 0.5.0

        """
        cdef int                      i
        cdef int                      j
        cdef CGappedSequence*         aligned
        cdef Sequence                 sequence
        cdef CSequence                cseq
        cdef vector[CSequence]        seqvec
        cdef vector[CGappedSequence*] gapvec
        cdef const float**            matrix
        cdef Alignment                alignment = Alignment.__new__(Alignment)
        cdef CFAMSA*                  famsa     = new CFAMSA(self._params)

        try:
            with nogil:
                famsa.adjustParams(profile1._msa.size() + profile2._msa.size())
                # copy score matrix weights
                self._copy_matrix(famsa)
                # align profiles
                if not famsa.alignProfiles(profile1._msa, profile2._msa):
                    raise RuntimeError("failed to align profiles")
                # take ownership of the final alignment
                for aligned in famsa.final_profile.data:
                    alignment._msa.emplace_back(aligned)
                famsa.final_profile.data.clear()
        finally:
            del famsa

        return alignment

    cpdef GuideTree build_tree(self, object sequences):
        """Build a tree from the given sequences.

        Arguments:
            sequences (iterable of `~pyfamsa.Sequence`): An iterable
                yielding the digitized sequences to build a tree from.

        Returns:
            `~pyfamsa.GuideTree`: The guide tree obtained from the sequences.

        """
        cdef size_t                            i
        cdef Sequence                          sequence
        cdef CSequence                         cseq
        cdef vector[CSequence]                 seqvec
        cdef vector[CSequence*]                ptrvec
        cdef shared_ptr[AbstractTreeGenerator] gen
        cdef vector[int]                       og2map
        cdef const float**                     matrix
        cdef CFAMSA*                           famsa    = new CFAMSA(self._params)
        cdef GuideTree                         tree     = GuideTree.__new__(GuideTree)

        try:
            # copy score matrix weights
            with nogil:
                self._copy_matrix(famsa)
            # copy the aligner input and record original order
            for i, sequence in enumerate(sequences):
                cseq = CSequence(sequence._cseq.get()[0])
                cseq.sequence_no = cseq.original_no = i
                seqvec.push_back(move(cseq))
            # sort sequences and record pointers
            if seqvec.size() > 0:
                famsa.sortAndExtendSequences(seqvec)
            for i in range(seqvec.size()):
                og2map.push_back(i)
                ptrvec.push_back(&seqvec.data()[i])
                tree._names.push_back(move(CSequence(seqvec[i].id, string(), i, NULL)))
            # remove duplicates and record sequence order
            if not self._params.keepDuplicates:
                famsa.removeDuplicates(ptrvec, og2map)
            for i in range(ptrvec.size()):
                ptrvec[i].sequence_no = i
            # generate tree
            if ptrvec.size() > 1:
                with nogil:
                    gen = famsa.createTreeGenerator(self._params)
                    gen.get().call(ptrvec, tree._tree.raw() )
            elif ptrvec.size() == 1:
                tree._tree.raw().push_back(node_t(-1, -1))
        finally:
            del famsa

        # return the result tree
        tree._tree.fromUnique(og2map)
        return tree


cdef class GuideTree:
    """A guide tree generated from several sequences.

    Guide trees are binary trees, stored as an array which contains the
    two children nodes for each node (or *-1* for leaves). Nodes can be
    iterated upon.

    """

    # --- Magic methods ------------------------------------------------------

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

    # --- Properties ---------------------------------------------------------

    @property
    def names(self):
        """`list` of `bytes`: The name of the leaves in the guide tree.
        """
        return [ <bytes> self._names[i].id for i in range(self._names.size()) ]

    # --- Methods ------------------------------------------------------------

    cpdef bytes dumps(self):
        """Dump the tree in Newick format into a `bytes` object.

        Returns:
            `bytes`: The tree in Newick format, stored as an ASCII string.

        """
        cdef string       out
        cdef NewickParser nw

        with nogil:
            nw.store(self._names, self._tree.raw(), out)

        return <bytes> out

    cpdef ssize_t dump(self, object file) except -1:
        """Write the tree into a file in Newick format.

        Arguments:
            file (`str`, `bytes, `os.PathLike` or file-like object): The
                path to a file, or a file-like object open in binary mode.

        """
        cdef string       out
        cdef NewickParser nw
        cdef object       mem
        cdef bytes        path = None

        # accept either a path or a file-like object
        if isinstance(file, (os.PathLike, str)):
            path = os.fsencode(file)
        elif isinstance(file, bytes):
            path = file
        else:
            path = None

        # dump the tree to a C++ string
        with nogil:
            nw.store(self._names, self._tree.raw(), out)

        # use a memoryview to avoid copy from the C++ string before writing
        mem = PyMemoryView_FromMemory(<char*> out.c_str(), out.size(), PyBUF_READ )

        # write the result to the file and return the number of bytes written
        if path is None:
            return file.write(mem)
        else:
            with open(path, "wb") as handle:
                return handle.write(mem)
