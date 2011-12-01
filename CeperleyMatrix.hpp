#ifndef _CEPERLEY_MATRIX_HPP
#define _CEPERLEY_MATRIX_HPP

#include <iostream>
#include <cmath>

#include <Eigen/Dense>
#include <boost/assert.hpp>

/**
 * O(N^2) method for keeping track of a determinant when only one row (or
 * column) changes in a given step.  Also known as the
 * Sherman-Morrison-Woodbury formula.  This class acts as a finite state
 * machine; that is, its methods should be called in a specific order.  See
 * next_step.
 */
template <typename T>
class CeperleyMatrix
{
private:
    /**
     * An enum for storing which operation should be called next on the object
     */
    enum NextStep {
        INITIALIZE,
        UPDATE,
        FINISH_ROW_UPDATE,
        FINISH_COLUMN_UPDATE
    };

    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat, invmat;
    T detrat, det;
    unsigned int pending_index; // refers to a row or column index
    NextStep next_step;

public:
    /**
     * Constructor for initializing a CeperleyMatrix from a square Eigen::Matrix
     */
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

    /**
     * It's sometimes useful to have a default constructor, but such objects
     * are useless until they are copied from an existing CeperleyMatrix.
     */
    CeperleyMatrix (void)
        : next_step(INITIALIZE)
        {
        }

    /**
     * Call this function to swap two rows in the matrix.  As a result, the
     * determinant will change sign.
     *
     * @param r1 Index of one row to be swapped
     * @param r2 Index of the row to the swapped with
     */
    void swap_rows (unsigned int r1, unsigned int r2)
        {
            BOOST_ASSERT(next_step == UPDATE);
            BOOST_ASSERT(r1 < mat.rows());
            BOOST_ASSERT(r2 < mat.rows());
            BOOST_ASSERT(r1 != r2);

            mat.row(r1).swap(mat.row(r2));
            invmat.col(r1).swap(invmat.col(r2));

            det = -det;
            //detrat = -detrat;
        }

    /**
     * Update a row in the matrix by replacing it with the given vector.
     *
     * This takes O(N) time.
     *
     * If the row can be represented roughly as a linear combination of the
     * existing rows, the determinant should become zero (in theory, but not
     * practice), and there will likely be large error from then on.
     *
     * After this is called, the new determinant is available, but no other
     * operations can be called until finish_row_update() is called.  We don't
     * do these steps here because they take O(N^2) time and are irrelevant if
     * the matrix is immediately thrown away (e.g. if the Monte Carlo step is
     * rejected).
     *
     * @param r index of the row to be updated
     *
     * @param row vector containing the values with which the row should be
     * replaced
     *
     * @see finish_row_update()
     * @see update_column()
     */
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

    /**
     * Update a column in the matrix by replacing it with the given vector.
     *
     * @see update_row()
     * @see finish_column_update()
     */
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

    /**
     * Finish a row update.  Must be called after update_row().
     *
     * This takes O(N^2) time.
     *
     * @see update_row()
     * @see finish_column_update()
     */
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

    /**
     * Finish a column update.  Must be called after update_column().
     *
     * @see update_column()
     * @see finish_row_update()
     */
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

    /**
     * Recalculates the matrix inverse and determinant from scratch.
     */
    void refresh_state (void)
        {
            BOOST_ASSERT(next_step == UPDATE);
            Eigen::FullPivLU<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > fullpivlu_decomposition(mat);
            invmat = fullpivlu_decomposition.inverse();
            det = fullpivlu_decomposition.determinant();
        }

    /**
     * Returns the matrix.
     */
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & get_matrix (void) const
        {
            BOOST_ASSERT(next_step != INITIALIZE);
            return mat;
        }

    /**
     * Returns the inverse matrix.
     */
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & get_inverse (void) const
        {
            BOOST_ASSERT(next_step == UPDATE);
            return invmat;
        }

    /**
     * Returns the determinant of the matrix.
     *
     * The determinant is always pre-calculated, so this runs in O(1) time.
     */
    T get_determinant (void) const
        {
            BOOST_ASSERT(next_step != INITIALIZE);
            return det;
        }

    /**
     * Computes and returns the inverse matrix error.
     *
     * This multiplies the matrix by its supposed inverse, and compares it to
     * the identity matrix.
     *
     * @return a nonnegative number which is the sum of the absolute error.
     */
    double compute_inverse_matrix_error (void) const
        {
            BOOST_ASSERT(next_step == UPDATE);
            return (mat * invmat - Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(mat.rows(), mat.cols())).array().abs().sum();
        }

    /**
     * Calculates and returns the determinant error.
     *
     * @return a ratio of the absolute error over the newly-calculated
     * determinant
     */
    double compute_relative_determinant_error (void) const
        {
            BOOST_ASSERT(next_step == UPDATE);
            T d = mat.fullPivLu().determinant();
            return abs((d - det) / d);
        }

    /**
     * Number of rows in the matrix
     *
     * (Since the matrix is square, this will always be equal to the number of
     * columns)
     *
     * @see cols()
     */
    unsigned int rows (void) const
        {
            BOOST_ASSERT(next_step != INITIALIZE);
            return mat.rows();
        }

    /**
     * Number of columns in the matrix
     *
     * (Since the matrix is square, this will always be equal to the number of
     * rows)
     *
     * @see rows()
     */
    unsigned int cols (void) const
        {
            BOOST_ASSERT(next_step != INITIALIZE);
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
