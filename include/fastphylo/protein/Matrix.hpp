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
 * on t. The decomposition itself is done once, in the constructor -
 * callers that need to evaluate exp(Q*t) more than once (any number of
 * times: a handful of Newton iterations, hundreds of bootstrap
 * replicates, or an unbounded number of calls from a long-lived
 * caller like fastphylo-py bindings) construct exactly one MatrixExpm
 * per model and reuse it - see calculate_ml_dists()'s doc comment
 * (ProtDistCalc.hpp) for why this is a hard requirement, not just an
 * optimization: nothing in this class's public API lets a caller
 * accidentally redo the decomposition per call the way a "build me
 * a fresh one from model_type" convenience function would invite.
 *
 * Two constructors, matching two different situations: the
 * Matrix+equilibrium-distribution one uses
 * Eigen::SelfAdjointEigenSolver on the symmetric similarity transform
 * S = D^0.5 * Q * D^-0.5 (D = diag(eq)), not the general (possibly
 * complex-eigenvalue) EigenSolver directly on Q - valid because every
 * rate matrix this constructor is meant for is time-reversible
 * (eq[i]*Q(i,j) == eq[j]*Q(j,i) by construction - true for all 5 named
 * protein models, see ModelMatrix.cpp; -F/--model-file custom matrices
 * aren't implemented yet, so this can't be violated by a real caller
 * today, but whoever wires up -F should validate the property then,
 * not assume it silently holds for arbitrary user input). Q's own
 * eigenvectors/inverse are recovered from S's orthogonal eigenvectors
 * U as V = D^-0.5*U, V^-1 = U^T*D^0.5 - no explicit matrix inverse
 * needed. Measured 3.15-3.8x faster than the general EigenSolver +
 * explicit .inverse() the other constructor still uses
 * (benchmarks/bench_symmetric_decomp.cpp), exact to machine precision
 * for the 4 exactly-reversible models (WAG/JTT/Dayhoff/LG) and to
 * 6.2e-5 for MVR (which has a small, ~1.2% relative, genuine detailed-
 * balance violation in its published data at one cell pair - well
 * inside this project's established 5e-4 correctness tolerance). The
 * Matrix-only constructor (no eq) keeps the general EigenSolver path -
 * see its own doc comment for why it still exists.
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
 * A second, `float`-precision copy of the decomposition (eigenvectors/
 * eigenvectors_inv/eigenvalues/Q) is cached alongside the `double` one
 * - likelihood_slope_curv()'s per-iteration arithmetic runs in `float`
 * (planning/fastprot_ml_speedup_implementation_plan.md: ~2x further
 * win on top of `double`-Eigen, error 3-5 orders of magnitude inside
 * the solver's convergence tolerance, measured on real data across
 * all 5 models). The decomposition itself stays `double`-only - it's
 * a one-time cost (not worth floating for speed) and the more
 * numerically sensitive step (eigendecomposition feeding a matrix
 * exponential); only the cheap, per-t evaluation is done in `float`.
 *
 * Also caches Q*Q (Q2_eigen_f(), `float` only - only the `float` hot
 * path uses it), computed once here rather than derived as
 * P''(t) = P'(t)*Q on every iteration: P'(t) and P''(t) can then both
 * be computed directly from P(t) (P'(t)=P(t)*Q, P''(t)=P(t)*Q²)
 * instead of P''(t) depending on P'(t) first - same total
 * multiplication count, but no dependency chain between the two.
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
    //! General-purpose decomposition (Eigen's general EigenSolver,
    //! discarding any imaginary part - matching this constructor's
    //! pre-2026-08-11 behavior exactly) - kept only for Matrix::expm()
    //! (the ED path, ExpectedDistance.cpp), which calls it on a Matrix
    //! with no equilibrium distribution available at that call site
    //! and isn't performance-critical (see
    //! planning/fastprot_ml_speedup_investigation_plan.md's scope
    //! notes). New callers that have eq available - which today means
    //! every named protein model - should use the constructor below
    //! instead: it's 3.15-3.8x faster and this one will silently give
    //! a wrong answer for a non-reversible Q where that one would
    //! correctly reject it (moot today since nothing not already known
    //! to be reversible reaches either constructor - see the class
    //! comment).
    explicit MatrixExpm(const Matrix &Q, double time_unit_scale = 1.0);
    //! Decomposes Q*time_unit_scale via the symmetrized approach
    //! described in the class comment - eq is Q's equilibrium
    //! distribution (needed for the D=diag(eq) similarity transform,
    //! not just a convenience). Q must be 20x20, eq must have 20
    //! entries (see class comment). time_unit_scale rescales Q before
    //! decomposing (default 1, i.e. no rescaling) - lets a caller
    //! change what one unit of the `t` it later passes to
    //! at()/at_eigen()/at_eigen_f() means, without needing a separate
    //! rescaling step at every call site that uses the result (see
    //! ProtDistCalc.cpp's build_ml_decomposition() for why: the
    //! published rate matrices are calibrated in a model-dependent
    //! scale, not the substitutions-per-site units fastprot's output
    //! uses).
    explicit MatrixExpm(const Matrix &Q, const DblVec &eq, double time_unit_scale = 1.0);
    //! Evaluates exp(Q*t) using the cached decomposition, as a Matrix
    //! - used by Matrix::expm() (the ED path, not the ML hot path
    //! this class was optimized for and not performance-critical -
    //! see planning/fastprot_ml_speedup_investigation_plan.md's scope
    //! notes on ExpectedDistance.cpp).
    Matrix at(double t) const;
    //! Evaluates exp(Q*t) using the cached decomposition, as a fixed-
    //! size Eigen matrix directly - not used by the ML hot path itself
    //! (see at_eigen_f()) but kept for the double-precision case in
    //! general (at()'s implementation, tests).
    Eigen::Matrix<double, 20, 20> at_eigen(double t) const;
    //! Q, cached in the same Eigen form used internally - avoids
    //! likelihood_slope_curv() re-converting it from Matrix every call.
    const Eigen::Matrix<double, 20, 20> &Q_eigen() const
    {
        return m_Q;
    }
    //! Evaluates exp(Q*t) using the cached float-precision decomposition
    //! - the fast path for likelihood_slope_curv()'s per-Newton-
    //! iteration hot loop.
    Eigen::Matrix<float, 20, 20> at_eigen_f(float t) const;
    //! Q, cached in float form.
    const Eigen::Matrix<float, 20, 20> &Q_eigen_f() const
    {
        return m_Q_f;
    }
    //! Q*Q, cached in float form - see class comment.
    const Eigen::Matrix<float, 20, 20> &Q2_eigen_f() const
    {
        return m_Q2_f;
    }

  private:
    Eigen::Matrix<double, 20, 20> m_eigenvectors;     // T
    Eigen::Matrix<double, 20, 20> m_eigenvectors_inv; // T^-1
    Eigen::Matrix<double, 20, 1> m_eigenvalues_real;  // real part of Q's eigenvalues
    Eigen::Matrix<double, 20, 20> m_Q;

    Eigen::Matrix<float, 20, 20> m_eigenvectors_f;
    Eigen::Matrix<float, 20, 20> m_eigenvectors_inv_f;
    Eigen::Matrix<float, 20, 1> m_eigenvalues_real_f;
    Eigen::Matrix<float, 20, 20> m_Q_f;
    Eigen::Matrix<float, 20, 20> m_Q2_f;
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
