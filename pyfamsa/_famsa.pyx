# coding: utf-8
# cython: language_level=3, linetrace=True

from libcpp.utility cimport move
from libcpp.vector cimport vector
from libcpp.string cimport string

from famsa.msa cimport CFAMSA
from famsa.core.io_service cimport IOService
from famsa.core.params cimport CParams
from famsa.core.sequence cimport CSequence, CGappedSequence
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

    def __init__(self):
        pass

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
