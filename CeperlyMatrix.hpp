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
	FINISH_ROW_UPDATE
    };

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat, invmat;
    T detrat, det;
    unsigned int pending_index;
    NextStep next_step;

public:
    CeperlyMatrix (const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &initial_mat)
	: mat(initial_mat),
	  invmat(mat.inverse()),
	  detrat(0),
	  det(mat.determinant()),
	  next_step(UPDATE_ROW)
	{
	    // if invmat.inverse() is not close to mat it probably means our
	    // orbitals are probably not linearly independent!
	    double inverse_error = compute_inverse_matrix_error();
	    if (inverse_error > .0001)
		std::cerr << "Warning: inverse matrix error of " << inverse_error << std::endl;
	}

    CeperlyMatrix (void)
	: next_step(INITIALIZE)
	{
	}

    void swap_rows (unsigned int r1, unsigned int r2)
	{
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    BOOST_ASSERT(r1 < mat.rows());
	    BOOST_ASSERT(r2 < mat.rows());
	    BOOST_ASSERT(r1 != r2);

	    mat.row(r1).swap(mat.row(r2));
	    invmat.col(r1).swap(invmat.col(r2));

	    det = -det;
	    detrat = -detrat;
	}

    void update_row (unsigned int r, const Eigen::Matrix<T, Eigen::Dynamic, 1> &row)
	{
	    BOOST_ASSERT(r < mat.rows());
	    BOOST_ASSERT(row.rows() == mat.cols());
	    BOOST_ASSERT(next_step == UPDATE_ROW);

	    // update matrix
	    mat.row(r) = row;
	    pending_index = r;

	    // calculate determinant ratio
	    detrat = mat.row(pending_index) * invmat.col(pending_index);
	    det *= detrat;

#ifdef CAREFUL
	    // check to make sure the row given doesn't already exist in the
	    // matrix, thus ruining things by setting the determinant to zero
	    for (unsigned int i = 0; i < mat.rows(); ++i) {
		if (i != pending_index) {
		    if (mat.row(i) == mat.row(pending_index))
			std::cerr << "!" << i << "," << pending_index << ' ' << detrat << std::endl;
		    BOOST_ASSERT(mat.row(i) != mat.row(pending_index));
		}
	    }
#endif

	    next_step = FINISH_ROW_UPDATE;
	}

    T get_determinant_ratio (void) const
	{
	    BOOST_ASSERT(next_step != INITIALIZE);
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

#ifdef CAREFUL
	    if (compute_inverse_matrix_error() > 1) {
		std::cerr << "Recomputing inverse due to large inverse matrix error of " << compute_inverse_matrix_error() << std::endl;
		refresh_state();
	    }
#endif
#ifdef CAREFUL
	    if (compute_relative_determinant_error() > .03)
		std::cerr << "large determinant error! " << compute_relative_determinant_error() << std::endl;
#endif
	}

    void refresh_state (void)
	{
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    invmat = mat.inverse();
	    det = mat.determinant();
	    // FIXME: adjust detrat accordingly
	}

    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & get_matrix (void) const
	{
	    BOOST_ASSERT(next_step != INITIALIZE);
	    return mat;
	}

    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & get_inverse (void) const
	{
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    return invmat;
	}

    T get_determinant (void) const
	{
	    BOOST_ASSERT(next_step != INITIALIZE);
	    return det;
	}

    double compute_inverse_matrix_error (void) const
	{
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    return (mat * invmat - Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(mat.rows(), mat.cols())).array().abs().sum();
	}

    double compute_relative_determinant_error (void) const
	{
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    T d = mat.determinant();
	    return abs((d - det) / d);
	}

    unsigned int rows (void) const
	{
	    return mat.rows();
	}

    unsigned int cols (void) const
	{
	    return rows();
	}
};

#endif
