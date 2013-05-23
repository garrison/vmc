#ifndef _VMC_TWO_BODY_JASTROW_FACTOR_HPP
#define _VMC_TWO_BODY_JASTROW_FACTOR_HPP

// we always want to include vmc-typedefs.hpp before including Eigen
#include "vmc-typedefs.hpp"
#include <Eigen/Dense>

#include "JastrowFactor.hpp"
#include "vmc-real-part.hpp"

/**
 * Two body Jastrow factor
 *
 * $J = \exp ( -\frac{1}{2} \sum_{ij} u_{ij} n_i n_j )$
 *
 * $n_i$ counts the number of particles regardless of species
 */
template <typename AmplitudeType>
class TwoBodyJastrowFactor : public JastrowFactor<AmplitudeType>
{
public:
    TwoBodyJastrowFactor (const Eigen::Matrix<typename RealPart<AmplitudeType>::type, Eigen::Dynamic, Eigen::Dynamic> &correlation_matrix);

private:
    virtual Big<AmplitudeType> compute_jastrow (const PositionArguments &r) const override;

    Eigen::Matrix<typename RealPart<AmplitudeType>::type, Eigen::Dynamic, Eigen::Dynamic> m_correlation_matrix;

public:
    // CYTHON-LIMITATION: http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html#c-left-values
    static inline void set_matrix_coeff (Eigen::Matrix<typename RealPart<AmplitudeType>::type, Eigen::Dynamic, Eigen::Dynamic> &mat,
                                         unsigned int row, unsigned int col, typename RealPart<AmplitudeType>::type value)
        {
            mat(row, col) = value;
        }
};

#define VMC_SUPPORTED_AMPLITUDE_TYPE(type) extern template class TwoBodyJastrowFactor<type>
#include "vmc-supported-types.hpp"
#undef VMC_SUPPORTED_AMPLITUDE_TYPE

#endif
