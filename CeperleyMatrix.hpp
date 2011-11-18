#ifndef _CEPERLEY_MATRIX_HPP
#define _CEPERLEY_MATRIX_HPP

#include <iostream>
#include <cmath>

#include <Eigen/Dense>
#include <boost/assert.hpp>

// O(N^2) method for keeping track of a determinant when only one row changes
// in a given step.  Also known as the Sherman-Morrison-Woodbury formula.  This
// class acts as a finite state machine.  See next_step.
template <typename T>
class CeperleyMatrix
{
private:
    enum NextStep {
        INITIALIZE,
        UPDATE,
        FINISH_ROW_UPDATE,
        FINISH_COLUMN_UPDATE
    };

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat, invmat;
    T detrat, det;
    unsigned int pending_index;
    NextStep next_step;

public:
    CeperleyMatrix (const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &initial_mat)
        : mat(initial_mat),
          detrat(0),
          next_step(UPDATE)
        {
            Eigen::FullPivLU<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > fullpivlu_decomposition(mat);
            invmat = fullpivlu_decomposition.inverse();
            det = fullpivlu_decomposition.determinant();

            // if there is significant inverse error it probably means our
            // orbitals are not linearly independent!
            double inverse_error = compute_inverse_matrix_error();
            if (inverse_error > .0001)
                std::cerr << "Warning: inverse matrix error of " << inverse_error << std::endl;
        }

    CeperleyMatrix (void)
        : next_step(INITIALIZE)
        {
        }

    void swap_rows (unsigned int r1, unsigned int r2)
        {
            BOOST_ASSERT(next_step == UPDATE);
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
            BOOST_ASSERT(next_step == UPDATE);

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

    void update_column (unsigned int c, const Eigen::Matrix<T, Eigen::Dynamic, 1> &col)
        {
            BOOST_ASSERT(c < mat.cols());
            BOOST_ASSERT(col.rows() == mat.rows());
            BOOST_ASSERT(next_step == UPDATE);

            // update matrix
            mat.col(c) = col;
            pending_index = c;

            // calculate determinant ratio
            detrat = invmat.row(pending_index) * mat.col(pending_index);
            det *= detrat;

#ifdef CAREFUL
            // check to make sure the column given doesn't already exist in the
            // matrix, thus ruining things by setting the determinant to zero
            for (unsigned int i = 0; i < mat.cols(); ++i) {
                if (i != pending_index) {
                    if (mat.col(i) == mat.col(pending_index))
                        std::cerr << "!" << i << "," << pending_index << ' ' << detrat << std::endl;
                    BOOST_ASSERT(mat.col(i) != mat.col(pending_index));
                }
            }
#endif

            next_step = FINISH_COLUMN_UPDATE;
        }

    T get_determinant_ratio (void) const
        {
            BOOST_ASSERT(next_step != INITIALIZE);
            return detrat;
        }

    void finish_row_update (void)
        {
            BOOST_ASSERT(next_step == FINISH_ROW_UPDATE);

            // implement equation (12) of Ceperley et al, correctly given as eqn (4.22)
            // of Kent's thesis http://www.ornl.gov/~pk7/thesis/thesis.ps.gz
            Eigen::Matrix<T, Eigen::Dynamic, 1> oldcol(invmat.col(pending_index));
            invmat -= ((invmat.col(pending_index) / detrat) * (mat.row(pending_index) * invmat)).eval();
            invmat.col(pending_index) = oldcol / detrat;

            next_step = UPDATE;

#ifdef CAREFUL
            be_careful();
#endif
        }

    void finish_column_update (void)
        {
            BOOST_ASSERT(next_step == FINISH_COLUMN_UPDATE);

            // same as above: update the inverse matrix
            Eigen::Matrix<T, Eigen::Dynamic, 1> oldrow(invmat.row(pending_index));
            invmat -= ((invmat * mat.col(pending_index)) * (invmat.row(pending_index) / detrat)).eval();
            invmat.row(pending_index) = oldrow / detrat;

            next_step = UPDATE;

#ifdef CAREFUL
            be_careful();
#endif
        }

    void refresh_state (void)
        {
            BOOST_ASSERT(next_step == UPDATE);
            Eigen::FullPivLU<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > fullpivlu_decomposition(mat);
            invmat = fullpivlu_decomposition.inverse();
            det = fullpivlu_decomposition.determinant();
            // FIXME: adjust detrat accordingly
        }

    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & get_matrix (void) const
        {
            BOOST_ASSERT(next_step != INITIALIZE);
            return mat;
        }

    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & get_inverse (void) const
        {
            BOOST_ASSERT(next_step == UPDATE);
            return invmat;
        }

    T get_determinant (void) const
        {
            BOOST_ASSERT(next_step != INITIALIZE);
            return det;
        }

    double compute_inverse_matrix_error (void) const
        {
            BOOST_ASSERT(next_step == UPDATE);
            return (mat * invmat - Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(mat.rows(), mat.cols())).array().abs().sum();
        }

    double compute_relative_determinant_error (void) const
        {
            BOOST_ASSERT(next_step == UPDATE);
            T d = mat.fullPivLu().determinant();
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

private:
#ifdef CAREFUL
    void be_careful (void)
        {
            if (compute_inverse_matrix_error() > 1) {
                std::cerr << "Recomputing inverse due to large inverse matrix error of " << compute_inverse_matrix_error() << std::endl;
                refresh_state();
            }

            if (compute_relative_determinant_error() > .03)
                std::cerr << "large determinant error! " << compute_relative_determinant_error() << std::endl;
        }
#endif
};

#endif
