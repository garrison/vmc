cdef extern from "<memory>" namespace "std":
    cdef cppclass auto_ptr[T]:
        auto_ptr()
        auto_ptr(T*)
        auto_ptr(auto_ptr&)

        T* get()
        T& operator*()
        #T* operator->()
        #void operator=()
        T* release()
        void reset(T*)
        void reset()
