from libc.stdint cimport int64_t

from famsa.tree cimport tree_structure


cdef extern from "tree/GuideTree.h" nogil:

    cdef cppclass GuideTree:
        GuideTree()

        tree_structure& raw()
        int getSequenceCount()
        int64_t calculateSackinIndex()
