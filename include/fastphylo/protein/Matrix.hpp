#pragma once

#include <stdexcept>
#include <vector>
#include <Eigen/Dense>

/*
 * Class for handling matrices.
 * Internal data is stored in column major format for easier use with LAPACK.
 * Element at matrix[1][2] is accessed with matrix(1,2);
 */
class Matrix;
using DblVec = std::vector<double>;
using MatVec = std::vector<Matrix>;

class Matrix
{
  public:
    //! Default constructor
    Matrix();
    //! Constructor for square matrix
    explicit Matrix(std::size_t size);
    //! Constructor for rectangular matrix
    Matrix(std::size_t rows, std::size_t cols);
    //! Creates full matrix from array
    Matrix(const double *, std::size_t, std::size_t);
    //! Creates diagonal matrix from std::vector
    Matrix(const DblVec &);
    //! Copy constructor
    Matrix(const Matrix &);
    //! Assignment operator
    Matrix &operator=(const Matrix &);
    //! Returns the row dimension
    std::size_t get_rows() const
    {
        return nr_rows;
    }
    //! Returns the column dimension
    std::size_t get_cols() const
    {
        return nr_cols;
    }
    //! Index operator (row, col)
    double &operator()(const std::size_t, const std::size_t);
    //! Const index operator (row, col)
    const double operator()(const std::size_t, const std::size_t) const;
    //! Index operator for raw data
    double &operator()(const std::size_t);
    //! Const index operator for raw data
    const double operator()(const std::size_t) const;
    //! Matrix multiplication, with possibility to transpose the matrices
    static Matrix mult(const Matrix &lhv, const Matrix &rhv, bool tr_left, bool tr_right);
    //! Matrix multiplication
    static Matrix mult(const Matrix &lhv, const Matrix &rhv);
    //! Multiplication with a diagonal matrix
    Matrix diag_mult(const DblVec &) const;
    //! Elementwise logarithm of the matrix
    void mlog();
    //! Calculates the matrix exponential
    Matrix expm() const;
    //! Multiplies the matrix with a scalar, and then calculates the matrix exponential
    MatVec expm(const DblVec &) const;
    //! Prints the matrix
    void printm() const;
    void printmi() const;
    //! Sums all the elements
    double sum() const;
    //! Sums the diagonal elements
    double sum_diag() const;

  private:
    //! Stores the matrix in column-major order
    DblVec m_data;
    //! The amount of rows
    std::size_t nr_rows;
    //! The columns
    std::size_t nr_cols;
};

/*!
 * A cached eigendecomposition of a square matrix Q, letting exp(Q*t) be
 * evaluated for many different t values without re-decomposing Q each
 * time: Q = T*diag(eigenvalues)*T^-1, so exp(Q*t) =
 * T*diag(exp(eigenvalues*t))*T^-1 - only this last, cheap step depends
 * on t. The decomposition itself (Eigen's EigenSolver) is done once,
 * in the constructor.
 *
 * Internally uses Eigen::Matrix<double,20,20> (fixed-size, stack-
 * allocated, no external BLAS dispatch) rather than the general
 * Matrix/LAPACK-dgemm_ path Matrix's other operations use - a 20x20
 * dense multiply through a general-purpose BLAS was found to be the
 * single largest cost in fastprot's maximum-likelihood pipeline
 * (planning/fastprot_ml_speedup_investigation_plan.md's Phase 1
 * profiling), and Eigen's fixed-size type avoids that dispatch
 * overhead entirely (planning/fastprot_ml_speedup_implementation_
 * plan.md). 20 is hardcoded (not the general N the rest of this file
 * supports) because every real caller passes a 20-amino-acid protein
 * rate matrix - confirmed by grep before this change, no other caller
 * exists; the constructor throws if given anything else.
 *
 * Also caches Q itself (Q_eigen()), in the same Eigen form, since
 * likelihood_slope_curv() (MaximumLikelihood.cpp) needs it on every
 * Newton-Raphson iteration for P'(t)/P''(t) and shouldn't have to
 * re-convert it from Matrix each call.
 *
 * Added because likelihood_slope_curv() (MaximumLikelihood.cpp)
 * evaluates exp(Q*t) at a new t on every Newton-Raphson iteration, for
 * the same Q, for every pair in a run - see phase0_audit.md's "ML
 * speedup round" for the profiling behind the once-per-run caching,
 * and the two planning docs above for the Eigen backend.
 */
class MatrixExpm
{
  public:
    //! Decomposes Q. Imaginary parts of Q's eigenvalues are discarded,
    //! matching expm()'s existing behavior (valid for the real,
    //! diagonalizable rate matrices this is used for). Q must be
    //! 20x20 (see class comment).
    explicit MatrixExpm(const Matrix &Q);
    //! Evaluates exp(Q*t) using the cached decomposition, as a Matrix
    //! - used by Matrix::expm() (the ED path, not the ML hot path
    //! this class was optimized for and not performance-critical -
    //! see planning/fastprot_ml_speedup_investigation_plan.md's scope
    //! notes on ExpectedDistance.cpp).
    Matrix at(double t) const;
    //! Evaluates exp(Q*t) using the cached decomposition, as a fixed-
    //! size Eigen matrix directly - the fast path for
    //! likelihood_slope_curv()'s per-Newton-iteration hot loop, no
    //! Matrix round-trip.
    Eigen::Matrix<double, 20, 20> at_eigen(double t) const;
    //! Q, cached in the same Eigen form used internally - avoids
    //! likelihood_slope_curv() re-converting it from Matrix every call.
    const Eigen::Matrix<double, 20, 20> &Q_eigen() const
    {
        return m_Q;
    }

  private:
    Eigen::Matrix<double, 20, 20> m_eigenvectors;     // T
    Eigen::Matrix<double, 20, 20> m_eigenvectors_inv; // T^-1
    Eigen::Matrix<double, 20, 1> m_eigenvalues_real;  // real part of Q's eigenvalues
    Eigen::Matrix<double, 20, 20> m_Q;
};

// Non-member functions
//! Elementwise matrix multiplication
Matrix elem_mult(const Matrix &lhv, const Matrix &rhv);
//! Elementwise matrix division
Matrix elem_div(const Matrix &lhv, const Matrix &rhv);

/*!
 *  Inlined index operator
 *  @param row The row
 *  @param col The column
 *  @return A reference to the element at position (row, col)
 */
inline double &Matrix::operator()(const std::size_t row, const std::size_t col)
{
    if (row >= get_rows() || col >= get_cols())
        throw std::out_of_range("Indexing outside of matrix");
    return m_data[col * nr_rows + row];
}

/*!
 *  Inlined const index operator
 *  @param row The row
 *  @param col The column
 *  @return A const copy of the element at position (row, col)
 */
inline const double Matrix::operator()(const std::size_t row, const std::size_t col) const
{
    if (row >= get_rows() || col >= get_cols())
        throw std::out_of_range("Indexing outside of matrix");
    return m_data[col * nr_rows + row];
}

/*!
 *  Inlined index operator
 *  @param ind Index in the representation of the matrix
 *  @return A reference to the element at position (ind) in the datastructure for the matrix
 */
inline double &Matrix::operator()(const std::size_t ind)
{
    if (ind >= get_cols() * get_rows())
        throw std::out_of_range("Indexing outside of matrix");
    return m_data[ind];
}

/*!
 *  Inlined const index operator
 *  @param ind Index in the representation of the matrix
 *  @return A const copy of the element at position (ind) in the datastructure for the matrix
 */
inline const double Matrix::operator()(const std::size_t ind) const
{
    if (ind >= get_cols() * get_rows())
        throw std::out_of_range("Indexing outside of matrix");
    return m_data[ind];
}
