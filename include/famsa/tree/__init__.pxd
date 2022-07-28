from libcpp.pair cimport pair
from libcpp.vector cimport vector

cdef extern from "tree/TreeDefs.h" nogil:

      ctypedef pair[int, int] node_t
      ctypedef vector[node_t] tree_structure

      cdef cppclass GT:
          pass

cdef extern from "tree/TreeDefs.h" namespace "GT" nogil:
    cdef enum Method:
        SLINK
        MST_Prim
        UPGMA
        UPGMA_modified
        NJ
        chained
        imported

    cdef enum Heuristic:
        None
        PartTree
        ClusterTree
