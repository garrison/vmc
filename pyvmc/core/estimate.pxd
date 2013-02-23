from libcpp.vector cimport vector
from pyvmc.core cimport complex_t

cdef extern from "RunningEstimate.hpp":
    cdef cppclass CppIntegerRunningEstimate "RunningEstimate<unsigned int>":
        double get_recent_result()
        double get_cumulative_result()
        unsigned int get_num_recent_values()
        unsigned int get_num_cumulative_values()

    cdef cppclass CppRealRunningEstimate "RunningEstimate<real_t>":
        double get_recent_result()
        double get_cumulative_result()
        unsigned int get_num_recent_values()
        unsigned int get_num_cumulative_values()

    cdef cppclass CppComplexRunningEstimate "RunningEstimate<amplitude_t>":
        complex_t get_recent_result()
        complex_t get_cumulative_result()
        unsigned int get_num_recent_values()
        unsigned int get_num_cumulative_values()

cdef extern from "BinnedEstimate.hpp":
    cdef cppclass CppIntegerBinnedSum "BinnedEstimate<unsigned int>::BinnedSum":
        double get_mean()
        double get_error()
        unsigned int get_num_bins()

    cdef cppclass CppRealBinnedSum "BinnedEstimate<real_t>::BinnedSum":
        double get_mean()
        double get_error()
        unsigned int get_num_bins()

    cdef cppclass CppComplexBinnedSum "BinnedEstimate<amplitude_t>::BinnedSum":
        complex_t get_mean()
        double get_error()
        unsigned int get_num_bins()

    cdef cppclass CppIntegerBinnedEstimate "BinnedEstimate<unsigned int>" (CppIntegerRunningEstimate):
        const vector[CppIntegerBinnedSum] & get_binlevel_data()

    cdef cppclass CppRealBinnedEstimate "BinnedEstimate<real_t>" (CppRealRunningEstimate):
        const vector[CppRealBinnedSum] & get_binlevel_data()

    cdef cppclass CppComplexBinnedEstimate "BinnedEstimate<amplitude_t>" (CppComplexRunningEstimate):
        const vector[CppComplexBinnedSum] & get_binlevel_data()

cdef extern from "BlockedEstimate.hpp":
    cdef cppclass CppIntegerBlockedEstimate "BlockedEstimate<unsigned int>" (CppIntegerBinnedEstimate):
        const vector[double] & get_block_averages()
        unsigned int get_measurements_per_block()

    cdef cppclass CppRealBlockedEstimate "BlockedEstimate<real_t>" (CppRealBinnedEstimate):
        const vector[double] & get_block_averages()
        unsigned int get_measurements_per_block()

    cdef cppclass CppComplexBlockedEstimate "BlockedEstimate<amplitude_t>" (CppComplexBinnedEstimate):
        const vector[complex_t] & get_block_averages()
        unsigned int get_measurements_per_block()

cdef Estimate_from_CppIntegerRunningEstimate(const CppIntegerRunningEstimate& cpp)
cdef Estimate_from_CppRealRunningEstimate(const CppRealRunningEstimate& cpp)
cdef Estimate_from_CppComplexRunningEstimate(const CppComplexRunningEstimate& cpp)

cdef Estimate_from_CppIntegerBinnedEstimate(const CppIntegerBinnedEstimate& cpp)
cdef Estimate_from_CppRealBinnedEstimate(const CppRealBinnedEstimate& cpp)
cdef Estimate_from_CppComplexBinnedEstimate(const CppComplexBinnedEstimate& cpp)

cdef Estimate_from_CppIntegerBlockedEstimate(const CppIntegerBlockedEstimate& cpp)
cdef Estimate_from_CppRealBlockedEstimate(const CppRealBlockedEstimate& cpp)
cdef Estimate_from_CppComplexBlockedEstimate(const CppComplexBlockedEstimate& cpp)
