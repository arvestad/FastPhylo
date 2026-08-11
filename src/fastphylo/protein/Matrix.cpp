#include "fastphylo/protein/Matrix.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <cstdio>

// LAPACK/BLAS functions - just dgemm_ now. MatrixExpm's decomposition
// (dgeev_/dgetrf_/dgetri_ previously) moved to Eigen's EigenSolver -
// see the MatrixExpm implementation below and
// planning/fastprot_ml_speedup_implementation_plan.md. Matrix::mult()
// (below) is unrelated to that change and keeps using dgemm_ - it has
// other, non-hot-path callers (ExpectedDistance.cpp, etc.) this plan
// deliberately didn't touch.
extern "C"
{
    /*!
     * SUBROUTINE DGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
     *   for info, see: http://www.netlib.org/blas/dgemm.f
     */
    void dgemm_(char *transa, char *transb, int *m, int *n, int *k, double *alpha, double *a, int *lda, double *b,
                int *ldb, double *beta, double *c, int *ldc);
}
/*!
 *  Default constructor, constructs a matrix with size (0,0)
 */
Matrix::Matrix() : nr_rows(0), nr_cols(0), m_data(0)
{
}
/*!
 *  Constructs a square matrix, filled with 0
 *  @param size The side of the square matrix
 */
Matrix::Matrix(std::size_t size) : nr_rows(size), nr_cols(size), m_data(size * size, 0)
{
}

/*!
 *  Constructs a matrix with size (rows, cols), filled with 0
 *  @param rows Amount of rows
 *  @param cols Amount of columns
 */
Matrix::Matrix(std::size_t rows, std::size_t cols) : nr_rows(rows), nr_cols(cols), m_data(cols * rows, 0)
{
}

/*!
 *  Constructs a matrix with size (rows, cols) and fills it with elements from
 *  an array. The array needs to be of size rows*cols
 *  @param array[] The elements in row-major order
 *  @param rows The amount of rows
 *  @param cols The amount of cols
 */
Matrix::Matrix(const double *array, std::size_t rows, std::size_t cols)
    : m_data(rows * cols, 0), nr_rows(rows), nr_cols(cols)
{
    for (int i = 0; i < cols; i++)
    {
        for (int j = 0; j < rows; j++)
        {
            m_data[(i * rows) + j] = array[(j * cols) + i];
        }
    }
}
/*! Constructs a diagonal matrix from a vector
 *  @param d A reference to a vector containing the diagonal elements
 */
Matrix::Matrix(const std::vector<double> &d) : nr_rows(d.size()), nr_cols(d.size()), m_data(d.size() * d.size())
{
    for (int i = 0; i < d.size(); i++)
    {
        m_data[(i * d.size()) + i] = d[i];
    }
}

/*! Creates a copy of a matrix
 * @param m The matrix to be copied
 */
Matrix::Matrix(const Matrix &m) : nr_rows(m.get_rows()), nr_cols(m.get_cols()), m_data(m.m_data)
{
}
/*! Assignment operator
 *  @param m The matrix to be
 *  */
Matrix &Matrix::operator=(const Matrix &m)
{
    if (&m != this)
    {
        nr_rows = m.get_rows();
        nr_cols = m.get_cols();
        m_data.clear();
        m_data.reserve(m.get_rows() * m.get_cols());
        for (int i = 0; i < m.get_rows() * m.get_cols(); i++)
        {
            m_data.push_back(m(i));
        }
    }
    return *this;
}
/*! Function that multiplies two matrices and returns a copy of the product
 *  @param lhv The left matrix
 *  @param rhv The right matrix
 *  @return a const copy with the product of lhv*rhv
 */
Matrix Matrix::mult(const Matrix &lhv, const Matrix &rhv)
{
    return mult(lhv, rhv, false, false);
}

/*!
 *  Static function that multiplies to matrices and possiby transposes them,
 *  returns a copy of the product
 *  @param lhv The left matrix
 *  @param rhv The right matrix
 *  @param tr_left true if the left matrix should be transposed, false otherwise
 *  @param tr_right true if the right matrix should be transposed, false otherwise
 *  @return A const copy with the product lhv*rhv
 */
Matrix Matrix::mult(const Matrix &lhv, const Matrix &rhv, bool tr_left, bool tr_right)
{
    if (lhv.get_cols() != rhv.get_rows())
    {
        throw std::out_of_range("Matrix dimensions doesn't agree");
    }

    char transl;
    char transr;
    int left_op_rows = lhv.get_rows();
    int right_op_cols = rhv.get_cols();
    int left_op_cols = lhv.get_cols();
    double dummy_one = 1;
    double dummy_zero = 0;
    // Matrix to store the result
    Matrix result(left_op_rows, right_op_cols);

    // T - the matrix should be transposed
    // N - the matrix should not be transposed
    tr_left ? transl = 'T' : transl = 'N';
    tr_right ? transr = 'T' : transr = 'N';

    // dgemm_ is a BLAS subroutine that multiplies two matrices
    // const_cast<double*> because dgemm_ will not change those arguments
    dgemm_(&transl, &transr, &left_op_rows, &right_op_cols, &left_op_cols, &dummy_one,
           const_cast<double *>(lhv.m_data.data()), &left_op_cols, const_cast<double *>(rhv.m_data.data()),
           &right_op_cols, &dummy_zero, result.m_data.data(), &left_op_rows);
    return result;
}

// Should this function be turned into a non-member too?
/*!
 *
 * @return A copy of the product
 */
Matrix Matrix::diag_mult(const std::vector<double> &diag) const
{
    Matrix temp(*this);

    for (int i = 0; i < get_cols(); i++)
    {
        for (int j = 0; j < get_rows(); j++)
        {
            temp(i, j) = (*this)((i * get_rows()) + j) * diag[i];
        }
    }

    return temp;
}
/*!
 * A function that applies ln() on every element in the matrix.
 */
void Matrix::mlog()
{
    int size = get_rows() * get_cols();
    for (int i = 0; i < size; i++)
    {
        (*this)(i) = log((*this)(i));
    }
}

/*!
 *  A const function that calculates the matrix exponential and returns a copy.
 *  The exponential is calculated by finding the eigenvalues and eigenvectors
 *  of the matrix, exponentiating the eigenvalues. The eigenvalues is stored in
 *  a matrix V, eigenvectors is stored in a matrix A, inv(A) is calculated.
 *  The product A*V*inv(A) is returned.
 *  @return A copy of the exponentiated matrix
 */
Matrix Matrix::expm() const
{
    return expm(DblVec(1, 1))[0];
}

/*!
 *  A const function that for every value in a vector calculates the matrix
 *  exponential of the matrix multiplied with that value
 *  The exponential is calculated by finding the eigenvalues and eigenvectors
 *  of the matrix, exponentiating the eigenvalues. The eigenvalues is stored in
 *  a matrix V, eigenvectors is stored in a matrix A, inv(A) is calculated.
 *  The product A*V*inv(A) is returned.
 *  @param s A vector with values to be multiplied with the matrix before
 *            the exponent is calculated.
 *  @return A vector with the exponential of the matrix multiplied with every
 *            value in s
 */
MatVec Matrix::expm(const DblVec &s) const
{
    MatrixExpm decomp(*this);
    MatVec result;
    result.reserve(s.size());
    for (double d : s)
    {
        result.push_back(decomp.at(d));
    }
    return result;
}

namespace
{
Eigen::Matrix<double, 20, 20> to_eigen20(const Matrix &m)
{
    Eigen::Matrix<double, 20, 20> e;
    for (int r = 0; r < 20; r++)
    {
        for (int c = 0; c < 20; c++)
        {
            e(r, c) = m(r, c);
        }
    }
    return e;
}
} // namespace

/*!
 *  Decomposes Q = T*diag(eigenvalues)*T^-1 once (the expensive part of
 *  expm()), so at()/at_eigen() can evaluate exp(Q*t) for many t values
 *  cheaply. See the class comment in Matrix.hpp.
 */
MatrixExpm::MatrixExpm(const Matrix &Q)
{
    if (Q.get_rows() != Q.get_cols())
    {
        throw std::out_of_range("Matrix needs to be square");
    }
    if (Q.get_rows() != 20)
    {
        // Eigen::Matrix<double,20,20>'s fixed size means this class no
        // longer supports arbitrary sizes the way the LAPACK-based
        // version did - confirmed by grep before this change that
        // every real caller passes a 20x20 protein rate matrix (see
        // the class comment in Matrix.hpp).
        throw std::out_of_range("MatrixExpm only supports 20x20 matrices (the protein rate matrices it was built "
                                 "for)");
    }

    m_Q = to_eigen20(Q);
    Eigen::EigenSolver<Eigen::Matrix<double, 20, 20>> solver(m_Q);
    m_eigenvectors = solver.eigenvectors().real();
    m_eigenvectors_inv = m_eigenvectors.inverse();
    m_eigenvalues_real = solver.eigenvalues().real();

    // float copy for likelihood_slope_curv()'s per-iteration hot path
    // - the decomposition itself stays double-only above (a one-time
    // cost, and the more numerically sensitive step); Q*Q is computed
    // here in double for accuracy before casting down, though it's a
    // one-time cost either way. See the class comment in Matrix.hpp.
    m_eigenvectors_f = m_eigenvectors.cast<float>();
    m_eigenvectors_inv_f = m_eigenvectors_inv.cast<float>();
    m_eigenvalues_real_f = m_eigenvalues_real.cast<float>();
    m_Q_f = m_Q.cast<float>();
    m_Q2_f = (m_Q * m_Q).cast<float>();
}

/*!
 *  Evaluates exp(Q*t) = T*diag(exp(eigenvalues*t))*T^-1 using the
 *  cached decomposition, as a fixed-size Eigen matrix - no
 *  decomposition, no heap allocation, no external BLAS dispatch.
 */
Eigen::Matrix<double, 20, 20> MatrixExpm::at_eigen(double t) const
{
    Eigen::Matrix<double, 20, 1> scale = (m_eigenvalues_real * t).array().exp();
    return m_eigenvectors * scale.asDiagonal() * m_eigenvectors_inv;
}

/*!
 *  Evaluates exp(Q*t) using the cached float-precision decomposition -
 *  the fast path for likelihood_slope_curv()'s per-Newton-iteration
 *  hot loop (MaximumLikelihood.cpp).
 */
Eigen::Matrix<float, 20, 20> MatrixExpm::at_eigen_f(float t) const
{
    Eigen::Matrix<float, 20, 1> scale = (m_eigenvalues_real_f * t).array().exp();
    return m_eigenvectors_f * scale.asDiagonal() * m_eigenvectors_inv_f;
}

/*!
 *  Evaluates exp(Q*t) as a Matrix - used by Matrix::expm() (the ED
 *  path via ExpectedDistance.cpp, not the ML hot path this class was
 *  optimized for and not performance-critical, see
 *  planning/fastprot_ml_speedup_investigation_plan.md's scope notes).
 */
Matrix MatrixExpm::at(double t) const
{
    Eigen::Matrix<double, 20, 20> result = at_eigen(t);
    Matrix m(20, 20);
    for (int r = 0; r < 20; r++)
    {
        for (int c = 0; c < 20; c++)
        {
            m(r, c) = result(r, c);
        }
    }
    return m;
}

/*!
 *  A function for printing the matrix
 */
void Matrix::printm() const
{
    for (int i = 0; i < get_rows(); i++)
    {
        for (int j = 0; j < get_cols(); j++)
        {
            printf(" %.20e ", (*this)(i, j));
        }
        std::cout << std::endl;
    }
}
void Matrix::printmi() const
{
    for (int i = 0; i < get_rows(); i++)
    {
        for (int j = 0; j < get_cols(); j++)
        {
            printf(" %.3lf ", (*this)(i, j));
        }
        std::cout << std::endl;
    }
}
/*!
 *  A const function that sums all the elements in the matrix.
 *  @return The sum of all the elements in the matrix
 */
double Matrix::sum() const
{
    return std::accumulate(m_data.begin(), m_data.end(), 0.0);
}

/*!
 *  A const function that sums all the elements on the diagonal of the matrix
 *  @return The sum of all the diagonal elements of the matrix
 */
double Matrix::sum_diag() const
{
    double sum = 0.0;
    if (get_rows() != get_cols())
    {
        throw std::out_of_range("Matrix needs to be square");
    }

    for (int i = 0; i < get_rows(); i++)
    {
        sum += (*this)(i, i);
    }

    return sum;
}
//--------------------------------------------------------------
// Non-member functions

/*!
 *  Non-member function that multiplies two matrices elementwise. The product is
 *  C(i,j) = A(i,j)*B(i,j). The matrices must have the same dimensions.
 *  @param lhv The matrix to the right
 *  @param rhv The matrix to the right
 *  @param A copy of the product
 */
Matrix elem_mult(const Matrix &lhv, const Matrix &rhv)
{
    if (lhv.get_rows() != rhv.get_rows() || lhv.get_cols() != rhv.get_cols())
    {
        throw std::out_of_range("Matrix dimensions must agree");
    }

    Matrix temp(lhv.get_rows(), rhv.get_cols());

    for (int i = 0; i < lhv.get_rows() * lhv.get_cols(); i++)
    {
        temp(i) = lhv(i) * rhv(i);
    }

    return temp;
}

/*!
 *  Non-member function that divides a matrix with another elementwise. A copy
 *  of the quotient C(i,j) = A(i,j)/B(i,j) is returned.
 *  @param lhv The dividend
 *  @param rhv The divisor
 *  @return The quotient
 */
Matrix elem_div(const Matrix &lhv, const Matrix &rhv)
{
    if (lhv.get_rows() != rhv.get_rows() || lhv.get_cols() != rhv.get_cols())
    {
        throw std::out_of_range("Matrix dimensions must agree");
    }

    Matrix temp(lhv.get_rows(), lhv.get_cols());

    for (int i = 0; i < lhv.get_rows() * lhv.get_cols(); i++)
    {
        temp(i) = lhv(i) / rhv(i);
    }

    return temp;
}
