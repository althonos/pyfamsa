from libcpp cimport bool

cdef extern from "utils/log.h" nogil:

    cdef cppclass Log:
        void enable()
        void disable()
        bool isEnabled()

        @staticmethod
        Log& getInstance(int level) except +


cdef extern from "utils/log.h" namespace "Log" nogil:
    const int LEVEL_NORMAL
    const int LEVEL_VERBOSE
    const int LEVEL_DEBUG
