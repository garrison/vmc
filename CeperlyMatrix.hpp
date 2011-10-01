#ifndef _CEPERLY_MATRIX_HPP
#define _CEPERLY_MATRIX_HPP

#include <Eigen/Dense>
#include <boost/assert.hpp>

// O(N^2) algorithm for finding ratio of new to previous determinant
// This class acts as a finite state machine.  See next_step.
class CeperlyMatrix
{
private:
    enum NextStep {
	UPDATE_ROW,
	CALCULATE_DETERMINANT_RATIO,
	FINISH_ROW_UPDATE
    };

    Eigen::MatrixXd mat, invmat;
    double detrat;
    int pending_index;
    NextStep next_step;

public:
    CeperlyMatrix (const Eigen::MatrixXd &initial_mat)
	: mat(initial_mat),
	  invmat(mat.inverse()),
	  next_step(UPDATE_ROW)
	{
	}

    void update_row (int r, const Eigen::VectorXd &row)
	{
	    BOOST_ASSERT(r >= 0 && r < mat.cols());
	    BOOST_ASSERT(row.rows() == mat.rows());
	    BOOST_ASSERT(next_step == UPDATE_ROW);

	    mat.row(r) = row;

	    pending_index = r;
	    next_step = CALCULATE_DETERMINANT_RATIO;
	}

    double calculate_determinant_ratio (void)
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
	    Eigen::VectorXd oldcol(invmat.col(pending_index));
	    invmat -= (invmat.col(pending_index) / detrat * mat.row(pending_index) * invmat).eval();
	    invmat.col(pending_index) = oldcol / detrat;

	    next_step = UPDATE_ROW;
	}

    void recompute_inverse (void)
	{
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    invmat = mat.inverse();
	}

    double compute_inverse_matrix_error (void) const
	{
	    // there is surely a more informative way to do this
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    Eigen::MatrixXd m(mat.inverse() - invmat);
	    return m.array().abs().sum();
	}

    Eigen::MatrixXd get (void) const
	{
	    return mat;
	}

    Eigen::MatrixXd get_inverse (void) const
	{
	    BOOST_ASSERT(next_step == UPDATE_ROW);
	    return invmat;
	}

private:
    // disable default constructor
    CeperlyMatrix (void);
};

#endif
