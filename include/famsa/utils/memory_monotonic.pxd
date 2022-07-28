from libcpp cimport bool

cdef extern from "utils/memory_monotonic.h" namespace "refresh" nogil:

    cdef cppclass memory_monotonic_base:
        memory_monotonic_base(size_t block_size, size_t alignment)


    cdef cppclass memory_monotonic_unsafe(memory_monotonic_base):
        memory_monotonic_unsafe()
        memory_monotonic_unsafe(size_t _block_size, size_t _alignment)

        bool deallocation_status() except +
        void* allocate(size_t size) except +
        void deallocate[T](T*& p) except +
        void freeze() except +
        void release() except +
        void release_freezed() except +


    cdef cppclass memory_monotonic_safe(memory_monotonic_base):
        memory_monotonic_safe()
        memory_monotonic_safe(size_t _block_size, size_t _alignment)

        bool deallocation_status() except +
        void* allocate(size_t size) except +
        void deallocate[T](T*& p) except +
        void freeze() except +
        void release() except +
        void release_freezed() except +
