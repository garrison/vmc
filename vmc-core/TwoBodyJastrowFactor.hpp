#ifndef _VMC_TWO_BODY_JASTROW_FACTOR_HPP
#define _VMC_TWO_BODY_JASTROW_FACTOR_HPP

#include <Eigen/Dense>

#include "JastrowFactor.hpp"

/**
 * Two body Jastrow factor
 *
 * $J = \exp ( -\frac{1}{2} \sum_{ij} u_{ij} n_i n_j )$
 *
 * $n_i$ counts the number of particles regardless of species
 */
class TwoBodyJastrowFactor : public JastrowFactor
{
public:
    TwoBodyJastrowFactor (const Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> &correlation_matrix);

private:
    Big<amplitude_t> compute_jastrow (const PositionArguments &r) const;

    Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> m_correlation_matrix;

public:
    // CYTHON-LIMITATION: http://docs.cython.org/src/userguide/wrapping_CPlusPlus.html#c-left-values
    static inline void set_matrix_coeff (Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic> &mat,
                                         unsigned int row, unsigned int col, real_t value)
        {
            mat(row, col) = value;
        }
};

#endif
