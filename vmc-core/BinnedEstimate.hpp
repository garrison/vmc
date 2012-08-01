#ifndef _BINNED_ESTIMATE_HPP
#define _BINNED_ESTIMATE_HPP

#include <vector>

#include "RunningEstimate.hpp"

template <typename T>
class BinnedEstimate : public RunningEstimate<T>
{
public:
    void add_value (T value)
        {
            RunningEstimate<T>::add_value(value);
            // fixme: perform the binning
        }
};

#endif
