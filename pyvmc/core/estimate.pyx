cdef complex_from_cpp(const complex_t& c):
    return complex(c.real(), c.imag())

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
    rv = Estimate_from_CppIntegerRunningEstimate(cpp)
    return rv

cdef Estimate_from_CppRealBinnedEstimate(const CppRealBinnedEstimate& cpp):
    rv = Estimate_from_CppRealRunningEstimate(cpp)
    return rv

cdef Estimate_from_CppComplexBinnedEstimate(const CppComplexBinnedEstimate& cpp):
    rv = Estimate_from_CppComplexRunningEstimate(cpp)
    return rv
