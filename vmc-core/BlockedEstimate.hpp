#ifndef _VMC_BLOCKED_ESTIMATE_HPP
#define _VMC_BLOCKED_ESTIMATE_HPP

#include <vector>

#include <boost/assert.hpp>

#include "BinnedEstimate.hpp"

/**
 * This estimate allows us to retrieve several different block averages, all at
 * a blocking level that is automatically adjusted over time.
 */
template <typename T>
class BlockedEstimate : public BinnedEstimate<T>
{
public:
    typedef typename BinnedEstimate<T>::BinnedSum BinnedSum;
    typedef typename BinnedEstimate<T>::result_t result_t;

    /**
     * Once max_blocks is reached, the existing blocks will be combined
     * pairwise into half as many blocks, thus doubling the number of
     * measurements per block.
     */
    BlockedEstimate (unsigned int max_blocks=128)
        : m_half_max_blocks(max_blocks / 2),
          m_measurements_per_block(1),
          m_outstanding_measurements_count(0),
          m_outstanding_sum(T(0))
        {
            BOOST_ASSERT(max_blocks >= 2);
        }

    void add_value (T value)
        {
            BinnedEstimate<T>::add_value(value);

            // if we have the maximum number of blocks already, re-pack them at
            // the next higher block level
            if (m_block_averages.size() == 2 * m_half_max_blocks) {
                for (unsigned int i = 0; i < m_half_max_blocks; ++i) {
                    m_block_averages[i] = (m_block_averages[2 * i] + m_block_averages[2 * i + 1]) / real_t(2);
                }
                m_block_averages.resize(m_half_max_blocks);
                m_measurements_per_block *= 2;
            }
            BOOST_ASSERT(m_block_averages.size() < 2 * m_half_max_blocks);

            // tally the measurement
            m_outstanding_sum += value;
            ++m_outstanding_measurements_count;

            // add another block if necessary
            if (m_outstanding_measurements_count == m_measurements_per_block) {
                m_block_averages.push_back(m_outstanding_sum / real_t(m_measurements_per_block));
                m_outstanding_sum = T(0);
                m_outstanding_measurements_count = 0;
            }
            BOOST_ASSERT(m_outstanding_measurements_count < m_measurements_per_block);
        }

    const std::vector<result_t> & get_block_averages (void) const
        {
            return m_block_averages;
        }

    unsigned int get_measurements_per_block (void) const
        {
            return m_measurements_per_block;
        }

private:
    unsigned int m_half_max_blocks;
    unsigned int m_measurements_per_block;
    unsigned int m_outstanding_measurements_count;
    T m_outstanding_sum;
    std::vector<result_t> m_block_averages;
};

#endif
