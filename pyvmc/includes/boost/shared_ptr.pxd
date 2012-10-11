cdef extern from "<boost/shared_ptr.hpp>" namespace "boost":
    cdef cppclass shared_ptr[T]:
        shared_ptr() nogil
        shared_ptr(T*) nogil
        void reset() nogil
        void reset(T*) nogil
        T& operator*() nogil
        T* get() nogil
