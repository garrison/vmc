#include <boost/assert.hpp>

#include "TwoBodyJastrowFactor.hpp"

// FIXME: if we can improve the JastrowFactor API a bit, we should be able to
// do this in O(N*m) time instead of O(N^2) time, where m is the number of
// particles moved.

template<typename AmplitudeType>
TwoBodyJastrowFactor<AmplitudeType>::TwoBodyJastrowFactor (const Eigen::Matrix<typename RealPart<AmplitudeType>::type, Eigen::Dynamic, Eigen::Dynamic> &correlation_matrix)
    : m_correlation_matrix(correlation_matrix)
{
    // ensure the matrix is square and symmetric
    BOOST_ASSERT(correlation_matrix.rows() == correlation_matrix.cols());
    BOOST_ASSERT(correlation_matrix == correlation_matrix.transpose());
}

template<typename AmplitudeType>
Big<AmplitudeType> TwoBodyJastrowFactor<AmplitudeType>::compute_jastrow (const PositionArguments &r) const
{
    BOOST_ASSERT(r.get_N_sites() == m_correlation_matrix.cols());

    typename RealPart<AmplitudeType>::type corr = 0;

    for (unsigned int species1 = 0; species1 < r.get_N_species(); ++species1) {
        for (unsigned int index1 = 0; index1 < r.get_N_filled(species1); ++index1) {
            const Particle particle1(index1, species1);
            unsigned int index2_start = index1 + 1;
            for (unsigned int species2 = species1; species2 < r.get_N_species(); ++species2) {
                for (unsigned int index2 = index2_start; index2 < r.get_N_filled(species2); ++index2) {
                    const Particle particle2(index2, species2);
                    corr += m_correlation_matrix(r[particle1], r[particle2]);
                }
                index2_start = 0;
            }
            // same-particle term, necessary for getting the on-site term correct
            corr += .5 * m_correlation_matrix(r[particle1], r[particle1]);
        }
    }

    // returns exp(-corr)
    return Big<AmplitudeType>(1, -corr);
}

#define VMC_SUPPORTED_AMPLITUDE_TYPE(type) template class TwoBodyJastrowFactor<type>
#include "vmc-supported-types.hpp"
