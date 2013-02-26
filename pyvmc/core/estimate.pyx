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

class RestoredBinnedSum(BinnedSum):
    pass

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

    def to_hdf5(self, *args, **kwargs):
        return _save_estimate_to_hdf5(self, *args, **kwargs)

    @staticmethod
    def from_hdf5(*args, **kwargs):
        return RestoredEstimate(*args, **kwargs)

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

class RestoredEstimate(Estimate):
    def __init__(self, estimate_group):
        # from RunningEstimate
        self.result = estimate_group["result"][...]
        self.num_values = estimate_group.attrs["num_measurements"]

        # from BinnedEstimate
        self.binlevel_data = [RestoredBinnedSum(*args)
                              for args in zip(estimate_group["binlevel_mean_data"],
                                              estimate_group["binlevel_error_data"],
                                              estimate_group["binlevel_nbins_data"])]

        # from BlockedEstimate
        self.block_averages = estimate_group["block_averages"][...]
        self.measurements_per_block = estimate_group.attrs["measurements_per_block"]

def _save_estimate_to_hdf5(estimate, estimate_group):
    # from RunningEstimate
    estimate_group.create_dataset("result", data=estimate.result)
    estimate_group.attrs["num_measurements"] = estimate.num_values

    # from BinnedEstimate
    estimate_group.create_dataset("binlevel_mean_data", data=[d.mean for d in estimate.binlevel_data])
    estimate_group.create_dataset("binlevel_error_data", data=[d.error for d in estimate.binlevel_data])
    estimate_group.create_dataset("binlevel_nbins_data", data=[d.nbins for d in estimate.binlevel_data])

    # from BlockedEstimate
    estimate_group.create_dataset("block_averages", data=estimate.block_averages)
    estimate_group.attrs["measurements_per_block"] = estimate.measurements_per_block
