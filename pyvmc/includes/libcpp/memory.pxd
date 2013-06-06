cdef extern from "<memory>" namespace "std":
    cdef cppclass unique_ptr[T]:
        unique_ptr() nogil
        unique_ptr(T*) nogil
        unique_ptr(unique_ptr&) nogil

        T* get() nogil
        T& operator*() nogil
        #T* operator->()
        T* release() nogil
        void reset(T*) nogil
        void reset() nogil

    cdef cppclass shared_ptr[T]:
        shared_ptr() nogil
        shared_ptr(T*) nogil
        void reset() nogil
        void reset(T*) nogil
        T& operator*() nogil
        T* get() nogil
