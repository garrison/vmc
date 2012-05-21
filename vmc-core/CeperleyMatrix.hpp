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

    /**
     * As long as the determinant remains below this value, the O(N^2) update
     * algorithm will be used.  However, if the determinant falls below this
     * value we will recalculate the inverse from scratch to fight numerical
     * error.  This also allows us to determine when the matrix is singular.
     */
    static const T ceperley_determinant_cutoff;

public:
    /**
     * Constructor for initializing a CeperleyMatrix from a square Eigen::Matrix
     */
    CeperleyMatrix (const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &initial_mat)
        : mat(initial_mat),
          next_step(UPDATE)
        {
            BOOST_ASSERT(initial_mat.rows() == initial_mat.cols());

            calculate_inverse();
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
     * This takes O(N) time, assuming the matrix did not become singular as a
     * result of the most recent update.
     *
     * After this is called, the new determinant is available, but no other
     * operations can be called until finish_row_update() is called.  We don't
     * update the inverse matrix here because it takes O(N^2) time and is
     * irrelevant if the matrix is immediately thrown away (e.g. if the Monte
     * Carlo step is rejected).
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

            // We calculate the new determinant using the
            // Sherman-Morrison-Woodbury formula.
            detrat = mat.row(pending_index) * invmat.col(pending_index);
            det *= detrat;

            // If the determinant has become sufficiently small, we assume the
            // matrix has become singular
            if (std::abs(det) < std::abs(ceperley_determinant_cutoff))
                det = T(0);

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

            // We calculate the new determinant using the
            // Sherman-Morrison-Woodbury formula.
            detrat = invmat.row(pending_index) * mat.col(pending_index);
            det *= detrat;

            // If the determinant has become sufficiently small, we assume the
            // matrix has become singular
            if (std::abs(det) < std::abs(ceperley_determinant_cutoff))
                det = T(0);

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

            BOOST_ASSERT(det != T(0));

            // implement equation (12) of Ceperley et al, correctly given
            // as eqn (4.22) of Kent's thesis
            // http://www.ornl.gov/~pk7/thesis/thesis.ps.gz
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

            BOOST_ASSERT(det != T(0));

            // same as above in finish_row_update(): update the inverse
            // matrix
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
            calculate_inverse();
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
            Eigen::FullPivLU<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > lu_decomposition(mat);
            if (lu_decomposition.isInvertible()) {
                T d = lu_decomposition.determinant();
                return std::abs((d - det) / d);
            } else {
                return std::abs(det);
            }
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
    void calculate_inverse (void)
        {
            Eigen::FullPivLU<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > lu_decomposition(mat);

            // fixme (?)
            lu_decomposition.setThreshold(lu_decomposition.threshold() * 10);

            if (!lu_decomposition.isInvertible()) {
                // oddly enough, lu_decomposition.determinant() is not
                // guaranteed to return zero, so we handle this case
                // explicitly.
                det = T(0);
            } else {
                det = lu_decomposition.determinant();
                invmat = lu_decomposition.inverse();

                // if there is significant inverse error it probably means our
                // orbitals are not linearly independent!
                double inverse_error = compute_inverse_matrix_error();
                if (inverse_error > .0001)
                    std::cerr << "Warning: inverse matrix error of " << inverse_error << std::endl;
            }
        }

#ifdef CAREFUL
    void be_careful (void) const
        {
            if (det != T(0) && compute_inverse_matrix_error() > 1)
                std::cerr << "Large inverse matrix error of " << compute_inverse_matrix_error() << std::endl;

            if (compute_relative_determinant_error() > .03)
                std::cerr << "large determinant error! " << compute_relative_determinant_error() << std::endl;
        }
#endif
};

#endif
