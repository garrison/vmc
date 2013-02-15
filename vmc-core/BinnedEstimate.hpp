#ifndef _BINNED_ESTIMATE_HPP
#define _BINNED_ESTIMATE_HPP

#include <vector>
#include <cmath>

#ifndef BOOST_NUMERIC_FUNCTIONAL_STD_COMPLEX_SUPPORT
#define BOOST_NUMERIC_FUNCTIONAL_STD_COMPLEX_SUPPORT
#endif

#include <boost/assert.hpp>
#include <boost/accumulators/statistics/moment.hpp>

#include "RunningEstimate.hpp"

static inline bool is_just_below_a_power_of_two (unsigned int x)
{
    return !(x & (x + 1));
}

template <typename T>
class BinnedEstimate : public RunningEstimate<T>
{
public:
    class BinnedSum
    {
    public:
        typedef typename RunningEstimate<T>::result_t result_t;

        /**
         * Returns the mean of all measurements considered at this binning level.
         */
        result_t get_mean (void) const
            {
                return boost::accumulators::mean(acc);
            }

        /**
         * Returns the statistical error at this binning level
         */
        result_t get_error (void) const
            {
                BOOST_ASSERT(boost::accumulators::count(acc) >= 2);
                const result_t mean = boost::accumulators::mean(acc);
                const result_t variance = boost::accumulators::moment<2>(acc) - (mean * mean);
                return std::sqrt(std::abs(variance) / result_t(boost::accumulators::count(acc) - 1));
            }

        /**
         * Returns the number of bins at this binning level
         */
        unsigned int get_num_bins (void) const
            {
                return boost::accumulators::count(acc);
            }

    private:
        BinnedSum (T current_sum_)
            : current_sum(current_sum_)
            {
            }

        typedef boost::accumulators::accumulator_set<T, boost::accumulators::stats<boost::accumulators::tag::mean, boost::accumulators::tag::moment<2>, boost::accumulators::tag::count> > accumulator_t;

        T current_sum;
        accumulator_t acc;

        friend class BinnedEstimate<T>;
    };

    void add_value (T value)
        {
            // create a new bin-level if necessary
            if (is_just_below_a_power_of_two(this->get_num_cumulative_values()))
                binlevel_data.push_back(BinnedSum(this->get_cumulative_total_value()));

            RunningEstimate<T>::add_value(value);

            // perform the binning
            BOOST_ASSERT(1u << binlevel_data.size() > this->get_num_cumulative_values());
            BOOST_ASSERT(1u << (binlevel_data.size() - 1) <= this->get_num_cumulative_values());

            for (unsigned int i = 0; i < binlevel_data.size(); ++i)
                binlevel_data[i].current_sum += value;

            for (unsigned int i = 0; ; ++i) {
                BOOST_ASSERT(i < binlevel_data.size());
                binlevel_data[i].acc(binlevel_data[i].current_sum / double(1 << i));
                binlevel_data[i].current_sum = 0;
                if (this->get_num_cumulative_values() & (1 << i))
                    break;
            }
        }

    const std::vector<BinnedSum> & get_binlevel_data (void) const
        {
            return binlevel_data;
        }

protected:
    // binlevel_data[n] puts 2^n consecutive measurements in the same bin
    std::vector<BinnedSum> binlevel_data;
};

#endif
