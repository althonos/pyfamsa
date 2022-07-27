from libcpp cimport bool
from libcpp.string cimport string
from libcpp.vector cimport vector

from famsa.core.sequence cimport CSequence, CGappedSequence
from famsa.utils.memory_monotonic cimport memory_monotonic_safe

cdef extern from "core/io_service.h" nogil:

    cdef cppclass IOService:
        @staticmethod
        size_t loadFasta(const string& filename, vector[CSequence]& sequences) except +
        # @staticmethod
        # size_t loadFasta(const string& filename, vector[CSequence]& sequences, memory_monotonic_safe* mma) except +
        @staticmethod
        bool saveAlignment(const string& filename, vector[CGappedSequence*]& sequences, int no_threads, int gzip_level) except +
