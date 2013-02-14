cdef complex_from_cpp(const complex_t& c):
    return complex(c.real(), c.imag())

class BinnedSum(object):
    def __init__(self, mean, variance):
        self.mean = mean
        self.variance = variance

class Estimate(object):
    result = None
    num_values = None
    recent_result = None
    num_recent_values = None
    binlevel_data = None

    def __str__(self):
        return "<Estimate: {}>".format(str(self.result))

cdef Estimate_from_CppIntegerRunningEstimate(const CppIntegerRunningEstimate& cpp):
    rv = Estimate()
    rv.result = cpp.get_cumulative_result()
    rv.num_values = cpp.get_num_cumulative_values()
    rv.recent_result = cpp.get_recent_result()
    rv.num_recent_values = cpp.get_num_recent_values()
    return rv

cdef Estimate_from_CppRealRunningEstimate(const CppRealRunningEstimate& cpp):
    rv = Estimate()
    rv.result = cpp.get_cumulative_result()
    rv.num_values = cpp.get_num_cumulative_values()
    rv.recent_result = cpp.get_recent_result()
    rv.num_recent_values = cpp.get_num_recent_values()
    return rv

cdef Estimate_from_CppComplexRunningEstimate(const CppComplexRunningEstimate& cpp):
    rv = Estimate()
    rv.result = complex_from_cpp(cpp.get_cumulative_result())
    rv.num_values = cpp.get_num_cumulative_values()
    rv.recent_result = complex_from_cpp(cpp.get_recent_result())
    rv.num_recent_values = cpp.get_num_recent_values()
    return rv

cdef Estimate_from_CppIntegerBinnedEstimate(const CppIntegerBinnedEstimate& cpp):
    cdef unsigned int i
    rv = Estimate_from_CppIntegerRunningEstimate(cpp)
    rv.binlevel_data = [BinnedSum(cpp.get_binlevel_data()[i].get_mean(),
                                  cpp.get_binlevel_data()[i].get_variance())
                        for i in range(cpp.get_binlevel_data().size() - 1)]
    return rv

cdef Estimate_from_CppRealBinnedEstimate(const CppRealBinnedEstimate& cpp):
    cdef unsigned int i
    rv = Estimate_from_CppRealRunningEstimate(cpp)
    rv.binlevel_data = [BinnedSum(cpp.get_binlevel_data()[i].get_mean(),
                                  cpp.get_binlevel_data()[i].get_variance())
                        for i in range(cpp.get_binlevel_data().size() - 1)]
    return rv

cdef Estimate_from_CppComplexBinnedEstimate(const CppComplexBinnedEstimate& cpp):
    cdef unsigned int i
    rv = Estimate_from_CppComplexRunningEstimate(cpp)
    rv.binlevel_data = [BinnedSum(complex_from_cpp(cpp.get_binlevel_data()[i].get_mean()),
                                  complex_from_cpp(cpp.get_binlevel_data()[i].get_variance()))
                        for i in range(cpp.get_binlevel_data().size() - 1)]
    return rv
