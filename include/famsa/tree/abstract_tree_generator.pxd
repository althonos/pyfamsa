from libcpp.vector cimport vector

from famsa.core cimport instruction_set_t
from famsa.core.sequence cimport CSequence
from famsa.tree cimport tree_structure


cdef extern from "tree/AbstractTreeGenerator.h" nogil:

    cdef cppclass AbstractTreeGenerator:
        AbstractTreeGenerator(int n_threads, instruction_set_t instruction_set) except +

        void call "operator()"(vector[CSequence*]& sequences, tree_structure& tree) except +
