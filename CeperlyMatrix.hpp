#ifndef _CEPERLY_MATRIX_HPP
#define _CEPERLY_MATRIX_HPP

#include <Eigen/Dense>
#include <boost/assert.hpp>

// O(N^2) algorithm for finding ratio of new to previous determinant
// This class acts as a finite state machine.  See next_step.
template <typename T>
class CeperlyMatrix
{
private:
    enum NextStep {
	UPDATE_ROW,
	CALCULATE_DETERMINANT_RATIO,
	FINISH_ROW_UPDATE
    };

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat, invmat;
    T detrat;
    int pending_index;
    NextStep next_step;

public:
    CeperlyMatrix (const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &initial_mat)
	: mat(initial_mat),
	  invmat(mat.inverse()),
	  next_step(UPDATE_ROW)
	{
	    // if invmat.inverse() is not close to mat it probably means our
	    // orbitals are probably not linearly independent!
	    BOOST_ASSERT((mat - invmat.inverse()).array().abs().sum() < .00000001);
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

    void recompute_inverse (void)
	{
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    invmat = mat.inverse();
	}

    T compute_inverse_matrix_error (void) const
	{
	    // there is surely a more informative way to do this
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> m(mat.inverse() - invmat);
	    return m.array().abs().sum();
	}

    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & get (void) const
	{
	    return mat;
	}

    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & get_inverse (void) const
	{
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    return invmat;
	}

private:
    // disable default constructor
    CeperlyMatrix (void);
};

#endif
