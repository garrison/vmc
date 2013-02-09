from pyvmc.core cimport complex_t

cdef extern from "RunningEstimate.hpp":
    cdef cppclass CppIntegerRunningEstimate "RunningEstimate<unsigned int>":
        float get_recent_result()
        float get_cumulative_result()
        unsigned int get_num_recent_values()
        unsigned int get_num_cumulative_values()

    cdef cppclass CppRealRunningEstimate "RunningEstimate<real_t>":
        float get_recent_result()
        float get_cumulative_result()
        unsigned int get_num_recent_values()
        unsigned int get_num_cumulative_values()

    cdef cppclass CppComplexRunningEstimate "RunningEstimate<amplitude_t>":
        complex_t get_recent_result()
        complex_t get_cumulative_result()
        unsigned int get_num_recent_values()
        unsigned int get_num_cumulative_values()

cdef extern from "BinnedEstimate.hpp":
    cdef cppclass CppIntegerBinnedEstimate "BinnedEstimate<unsigned int>" (CppIntegerRunningEstimate):
        pass

    cdef cppclass CppRealBinnedEstimate "BinnedEstimate<real_t>" (CppRealRunningEstimate):
        pass

    cdef cppclass CppComplexBinnedEstimate "BinnedEstimate<amplitude_t>" (CppComplexRunningEstimate):
        pass
