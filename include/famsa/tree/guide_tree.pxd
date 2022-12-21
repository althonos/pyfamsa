from libc.stdint cimport int64_t
from libcpp.vector cimport vector

from famsa.tree cimport tree_structure


cdef extern from "tree/GuideTree.h" nogil:

    cdef cppclass GuideTree:
        GuideTree()

        tree_structure& raw()
        int getSequenceCount()
        int64_t calculateSackinIndex()

        void toUnique(const vector[int]& original2unique, int n_uniques) except +
        void fromUnique(const vector[int]& original2unique) except +