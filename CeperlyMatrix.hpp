#ifndef _CEPERLY_MATRIX_HPP
#define _CEPERLY_MATRIX_HPP

#include <iostream>
#include <cmath>

#include <Eigen/Dense>
#include <boost/assert.hpp>

// O(N^2) algorithm for finding ratio of new to previous determinant
// This class acts as a finite state machine.  See next_step.
template <typename T>
class CeperlyMatrix
{
private:
    enum NextStep {
	INITIALIZE,
	UPDATE_ROW,
	CALCULATE_DETERMINANT_RATIO,
	FINISH_ROW_UPDATE
    };

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat, invmat;
    T detrat, det;
    int pending_index;
    NextStep next_step;

public:
    CeperlyMatrix (const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &initial_mat)
	: mat(initial_mat),
	  invmat(mat.inverse()),
	  det(mat.determinant()),
	  next_step(UPDATE_ROW)
	{
	    // if invmat.inverse() is not close to mat it probably means our
	    // orbitals are probably not linearly independent!
	    double inverse_error = compute_inverse_matrix_error();
	    if (inverse_error > .0000000001)
		std::cerr << "Warning: inverse matrix error of " << inverse_error << std::endl;
	}

    CeperlyMatrix (void)
	: next_step(INITIALIZE)
	{
	}

    void swap_rows (int r1, int r2)
	{
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    BOOST_ASSERT(r1 >= 0 && r1 < mat.rows());
	    BOOST_ASSERT(r2 >= 0 && r2 < mat.rows());
	    BOOST_ASSERT(r1 != r2);

	    mat.row(r1).swap(mat.row(r2));
	    invmat.col(r1).swap(invmat.col(r2));

	    det = -det;
	}

    void update_row (int r, const Eigen::Matrix<T, Eigen::Dynamic, 1> &row)
	{
	    BOOST_ASSERT(r >= 0 && r < mat.rows());
	    BOOST_ASSERT(row.rows() == mat.cols());
	    BOOST_ASSERT(next_step == UPDATE_ROW);

	    mat.row(r) = row;

	    pending_index = r;
	    next_step = CALCULATE_DETERMINANT_RATIO;
	}

    T calculate_determinant_ratio (void)
	{
	    BOOST_ASSERT(next_step == CALCULATE_DETERMINANT_RATIO);

	    detrat = mat.row(pending_index) * invmat.col(pending_index);
	    det *= detrat;
#if 0
	    std::cerr << det << ' ' << mat.determinant() << std::endl;
#endif

	    next_step = FINISH_ROW_UPDATE;
	    return detrat;
	}

    void finish_row_update (void)
	{
	    BOOST_ASSERT(next_step == FINISH_ROW_UPDATE);

	    // implement equation (12) of Ceperly et al, correctly given as eqn (4.22)
	    // of Kent's thesis http://www.ornl.gov/~pk7/thesis/thesis.ps.gz
	    Eigen::Matrix<T, Eigen::Dynamic, 1> oldcol(invmat.col(pending_index));
	    invmat -= (invmat.col(pending_index) / detrat * mat.row(pending_index) * invmat).eval();
	    invmat.col(pending_index) = oldcol / detrat;

	    next_step = UPDATE_ROW;
	}

    void refresh_state (void)
	{
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    invmat = mat.inverse();
	    det = mat.determinant();
	}

    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & get_matrix (void) const
	{
	    return mat;
	}

    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & get_inverse (void) const
	{
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    return invmat;
	}

    T get_determinant (void) const
	{
#ifdef DEBUG
	    static unsigned int z = 0;
	    if (++z % 541 == 0) // use some prime number here
		std::cerr << "deterr " << compute_relative_determinant_error() << " (" << abs(det) << ")" << std::endl;
#endif
	    return det;
	}

    double compute_inverse_matrix_error (void) const
	{
	    // there is surely a more informative way to do this
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    return (mat * invmat - Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(mat.rows(), mat.cols())).array().abs().sum();
	}

    double compute_relative_determinant_error (void) const
	{
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    T d = mat.determinant();
	    return abs((d - det) / d);
	}
};

#endif
