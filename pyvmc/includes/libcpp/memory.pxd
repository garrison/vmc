cdef extern from "<memory>" namespace "std":
    cdef cppclass auto_ptr[T]:
        auto_ptr() nogil
        auto_ptr(T*) nogil
        auto_ptr(auto_ptr&) nogil

        T* get() nogil
        T& operator*() nogil
        #T* operator->()
        #void operator=()
        T* release() nogil
        void reset(T*) nogil
        void reset() nogil
