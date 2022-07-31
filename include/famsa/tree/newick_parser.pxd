from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from famsa.core.sequence cimport CSequence
from famsa.tree cimport tree_structure


cdef extern from "tree/NewickParser.h" nogil:

    cdef cppclass NewickParser:
        NewickParser()

        void parse(
            const vector[CSequence]& sequences,
            const string& description,
            tree_structure& guideTree,
        ) except +

        void store(
            const vector[CSequence]& sequences,
            const tree_structure& guideTree,
            string& description,
        ) except +
