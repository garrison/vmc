cdef complex_from_cpp(const complex_t& c):
    if c.imag() == 0:
        return c.real()
    else:
        return complex(c.real(), c.imag())

class BinnedSum(object):
    def __init__(self, mean, error, nbins):
        self.mean = mean
        self.error = error
        self.nbins = nbins

class Estimate(object):
    # from RunningEstimate
    result = None
    num_values = None

    # from BinnedEstimate
    binlevel_data = None

    # from BlockedEstimate
    block_averages = None
    measurements_per_block = None

    def __str__(self):
        return "<Estimate: {}>".format(str(self.result))

cdef Estimate_from_CppIntegerRunningEstimate(const CppIntegerRunningEstimate& cpp):
    rv = Estimate()
    rv.result = cpp.get_cumulative_result()
    rv.num_values = cpp.get_num_cumulative_values()
    return rv

cdef Estimate_from_CppRealRunningEstimate(const CppRealRunningEstimate& cpp):
    rv = Estimate()
    rv.result = cpp.get_cumulative_result()
    rv.num_values = cpp.get_num_cumulative_values()
    return rv

cdef Estimate_from_CppComplexRunningEstimate(const CppComplexRunningEstimate& cpp):
    rv = Estimate()
    rv.result = complex_from_cpp(cpp.get_cumulative_result())
    rv.num_values = cpp.get_num_cumulative_values()
    return rv

cdef Estimate_from_CppIntegerBinnedEstimate(const CppIntegerBinnedEstimate& cpp):
    cdef unsigned int i
    rv = Estimate_from_CppIntegerRunningEstimate(cpp)
    rv.binlevel_data = [BinnedSum(cpp.get_binlevel_data()[i].get_mean(),
                                  cpp.get_binlevel_data()[i].get_error(),
                                  cpp.get_binlevel_data()[i].get_num_bins())
                        for i in range(cpp.get_binlevel_data().size() - 1)]
    return rv

cdef Estimate_from_CppRealBinnedEstimate(const CppRealBinnedEstimate& cpp):
    cdef unsigned int i
    rv = Estimate_from_CppRealRunningEstimate(cpp)
    rv.binlevel_data = [BinnedSum(cpp.get_binlevel_data()[i].get_mean(),
                                  cpp.get_binlevel_data()[i].get_error(),
                                  cpp.get_binlevel_data()[i].get_num_bins())
                        for i in range(cpp.get_binlevel_data().size() - 1)]
    return rv

cdef Estimate_from_CppComplexBinnedEstimate(const CppComplexBinnedEstimate& cpp):
    cdef unsigned int i
    rv = Estimate_from_CppComplexRunningEstimate(cpp)
    rv.binlevel_data = [BinnedSum(complex_from_cpp(cpp.get_binlevel_data()[i].get_mean()),
                                  cpp.get_binlevel_data()[i].get_error(),
                                  cpp.get_binlevel_data()[i].get_num_bins())
                        for i in range(cpp.get_binlevel_data().size() - 1)]
    return rv

cdef Estimate_from_CppIntegerBlockedEstimate(const CppIntegerBlockedEstimate& cpp):
    cdef unsigned int i
    rv = Estimate_from_CppIntegerBinnedEstimate(cpp)
    rv.block_averages = [cpp.get_block_averages()[i]
                         for i in range(cpp.get_block_averages().size())]
    rv.measurements_per_block = cpp.get_measurements_per_block()
    return rv

cdef Estimate_from_CppRealBlockedEstimate(const CppRealBlockedEstimate& cpp):
    cdef unsigned int i
    rv = Estimate_from_CppRealBinnedEstimate(cpp)
    rv.block_averages = [cpp.get_block_averages()[i]
                         for i in range(cpp.get_block_averages().size())]
    rv.measurements_per_block = cpp.get_measurements_per_block()
    return rv

cdef Estimate_from_CppComplexBlockedEstimate(const CppComplexBlockedEstimate& cpp):
    cdef unsigned int i
    rv = Estimate_from_CppComplexBinnedEstimate(cpp)
    rv.block_averages = [complex_from_cpp(cpp.get_block_averages()[i])
                         for i in range(cpp.get_block_averages().size())]
    rv.measurements_per_block = cpp.get_measurements_per_block()
    return rv
