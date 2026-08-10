// Micro-benchmark for Q1 (fastprot ML speedup investigation plan,
// planning/fastprot_ml_speedup_investigation_plan.md): current
// Matrix (heap-allocated, dgemm_-based) vs. Eigen fixed-size
// Matrix<double,20,20> / Matrix<float,20,20>, on the exact operation
// Phase 1's profiling found dominates the ML pipeline today - the
// per-Newton-iteration exp(Q*t) evaluation in MatrixExpm::at(), i.e.
// a 20x20 column-scale followed by a 20x20 dense multiply.
//
// This is a standalone, exploratory benchmark (not wired into the
// CMake build) - Eigen is not yet a project dependency; this is what
// decides whether it should become one. The LAPACK decomposition
// (dgeev_/dgetrf_/dgetri_) is duplicated here rather than reaching
// into MatrixExpm's private fields, so production headers don't need
// a benchmark-only accessor added just for this investigation - it
// mirrors MatrixExpm's constructor (Matrix.cpp) exactly.
//
// Build manually:
//   c++ -O2 -std=c++17 -I../include \
//     -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 \
//     bench_matrix_backend.cpp \
//     ../src/fastphylo/protein/Matrix.cpp \
//     ../src/fastphylo/protein/ModelMatrix.cpp \
//     -framework Accelerate -o bench_matrix_backend
//   ./bench_matrix_backend

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "fastphylo/protein/Matrix.hpp"
#include "fastphylo/protein/ModelMatrix.hpp"

extern "C"
{
    void dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi, double *vl, int *ldvl,
                double *vr, int *ldvr, double *work, int *lwork, int *info);
    void dgetri_(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
    void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
}

using namespace std::chrono;

namespace
{

constexpr std::size_t N = 20;

struct Stats
{
    double median_us;
    double q1_us;
    double q3_us;
};

Stats summarize(std::vector<double> &samples_us)
{
    std::sort(samples_us.begin(), samples_us.end());
    std::size_t n = samples_us.size();
    Stats s;
    s.median_us = samples_us[n / 2];
    s.q1_us = samples_us[n / 4];
    s.q3_us = samples_us[(3 * n) / 4];
    return s;
}

template <typename Fn> Stats time_calls(Fn fn, int warmup, int reps)
{
    for (int i = 0; i < warmup; i++)
    {
        fn();
    }
    std::vector<double> samples;
    samples.reserve(reps);
    for (int i = 0; i < reps; i++)
    {
        auto t0 = high_resolution_clock::now();
        fn();
        auto t1 = high_resolution_clock::now();
        samples.push_back(duration<double, std::micro>(t1 - t0).count());
    }
    return summarize(samples);
}

void report(const std::string &name, const Stats &old_s, const Stats &new_s)
{
    std::cout << name << "," << old_s.median_us << "," << old_s.q1_us << "," << old_s.q3_us << "," << new_s.median_us
               << "," << new_s.q1_us << "," << new_s.q3_us << "," << (old_s.median_us / new_s.median_us) << "\n";
}

// Mirrors MatrixExpm's constructor (Matrix.cpp) exactly - decomposes
// Q = T*diag(eigenvalues)*T^-1 once via LAPACK.
struct LapackDecomp
{
    Matrix eigenvectors;
    Matrix eigenvectors_inv;
    DblVec eigenvalues_real;
};

LapackDecomp lapack_decompose(const Matrix &Q)
{
    int size = static_cast<int>(Q.get_rows());
    LapackDecomp d;
    d.eigenvalues_real.assign(size, 0);
    DblVec eg_val_im(size, 0);
    std::array<double, 1> dummy{};
    int dummy_size = 1;
    std::array<int, 1> info{};
    char n = 'N';
    char v = 'V';
    std::array<double, 1> workspace_size{};
    int w_query = -1;

    DblVec data(size * size);
    for (int i = 0; i < size * size; i++)
    {
        data[i] = Q(i);
    }

    d.eigenvectors = Matrix(size, size);
    std::vector<double> eigvecs_raw(size * size, 0);
    dgeev_(&n, &v, &size, data.data(), &size, d.eigenvalues_real.data(), eg_val_im.data(), dummy.data(), &dummy_size,
           eigvecs_raw.data(), &size, workspace_size.data(), &w_query, info.data());
    DblVec workspace_vec(static_cast<int>(workspace_size[0]), 0);
    int w_size = static_cast<int>(workspace_size[0]);

    // Re-copy Q's data since the workspace-query call may have
    // touched `data` (LAPACK's contract only guarantees A is
    // unmodified on a pure workspace query, but re-copying is cheap
    // and removes any doubt).
    for (int i = 0; i < size * size; i++)
    {
        data[i] = Q(i);
    }
    dgeev_(&n, &v, &size, data.data(), &size, d.eigenvalues_real.data(), eg_val_im.data(), dummy.data(), &dummy_size,
           eigvecs_raw.data(), &size, workspace_vec.data(), &w_size, info.data());
    for (int i = 0; i < size * size; i++)
    {
        d.eigenvectors(i) = eigvecs_raw[i];
    }

    d.eigenvectors_inv = d.eigenvectors;
    std::vector<int> ipiv(size);
    std::vector<double> inv_raw(size * size);
    for (int i = 0; i < size * size; i++)
    {
        inv_raw[i] = d.eigenvectors(i);
    }
    dgetrf_(&size, &size, inv_raw.data(), &size, ipiv.data(), info.data());
    dgetri_(&size, inv_raw.data(), &size, ipiv.data(), workspace_size.data(), &w_query, info.data());
    DblVec workspace_vec2(static_cast<int>(workspace_size[0]));
    w_size = static_cast<int>(workspace_size[0]);
    dgetri_(&size, inv_raw.data(), &size, ipiv.data(), workspace_vec2.data(), &w_size, info.data());
    for (int i = 0; i < size * size; i++)
    {
        d.eigenvectors_inv(i) = inv_raw[i];
    }

    return d;
}

using EMatd = Eigen::Matrix<double, N, N>;
using EMatf = Eigen::Matrix<float, N, N>;
template <typename Scalar> using EMat = Eigen::Matrix<Scalar, N, N>;

EMatd to_eigen(const Matrix &m)
{
    EMatd e;
    for (std::size_t r = 0; r < N; r++)
    {
        for (std::size_t c = 0; c < N; c++)
        {
            e(static_cast<int>(r), static_cast<int>(c)) = m(r, c);
        }
    }
    return e;
}

Matrix to_matrix(const EMatd &e)
{
    Matrix m(N, N);
    for (std::size_t r = 0; r < N; r++)
    {
        for (std::size_t c = 0; c < N; c++)
        {
            m(r, c) = e(static_cast<int>(r), static_cast<int>(c));
        }
    }
    return m;
}

// Naive triple-loop 20x20x20 multiply - the simplest possible
// baseline, no library involved at all.
Matrix naive_mult(const Matrix &a, const Matrix &b)
{
    Matrix result(N, N);
    for (std::size_t i = 0; i < N; i++)
    {
        for (std::size_t j = 0; j < N; j++)
        {
            double sum = 0;
            for (std::size_t k = 0; k < N; k++)
            {
                sum += a(i, k) * b(k, j);
            }
            result(i, j) = sum;
        }
    }
    return result;
}

// Mirrors MatrixExpm::at() exactly (Matrix.cpp).
Matrix expm_at_current(const Matrix &eigenvectors, const Matrix &eigenvectors_inv, const DblVec &eigenvalues_real,
                        double t)
{
    Matrix left(N, N);
    for (std::size_t col = 0; col < N; col++)
    {
        double scale = exp(eigenvalues_real[col] * t);
        for (std::size_t row = 0; row < N; row++)
        {
            left(row, col) = eigenvectors(row, col) * scale;
        }
    }
    return Matrix::mult(left, eigenvectors_inv);
}

template <typename Scalar>
EMat<Scalar> expm_at_eigen(const EMat<Scalar> &eigenvectors, const EMat<Scalar> &eigenvectors_inv,
                            const Eigen::Matrix<Scalar, N, 1> &eigenvalues_real, Scalar t)
{
    Eigen::Matrix<Scalar, N, 1> scale = (eigenvalues_real * t).array().exp();
    EMat<Scalar> left = eigenvectors * scale.asDiagonal();
    return left * eigenvectors_inv;
}

double max_abs_diff(const Matrix &a, const Matrix &b)
{
    double worst = 0;
    for (std::size_t r = 0; r < N; r++)
    {
        for (std::size_t c = 0; c < N; c++)
        {
            worst = std::max(worst, std::fabs(a(r, c) - b(r, c)));
        }
    }
    return worst;
}

} // namespace

int main()
{ // NOLINT(bugprone-exception-escape)
    const int WARMUP = 3;
    const int REPS = 21;

    // Real rate matrix (WAG), not a synthetic one - the actual shape
    // of data likelihood_deriv()/MatrixExpm::at() operate on.
    Matrix Q = get_model_matrix(wag);
    LapackDecomp lap = lapack_decompose(Q);

    Eigen::EigenSolver<EMatd> es(to_eigen(Q));
    EMatd eig_vectors_d = es.eigenvectors().real();
    EMatd eig_vectors_inv_d = eig_vectors_d.inverse();
    Eigen::Matrix<double, N, 1> eig_values_d = es.eigenvalues().real();

    EMatf eig_vectors_f = eig_vectors_d.cast<float>();
    EMatf eig_vectors_inv_f = eig_vectors_inv_d.cast<float>();
    Eigen::Matrix<float, N, 1> eig_values_f = eig_values_d.cast<float>();

    const double t = 0.35; // representative mid-range branch length

    // --- Correctness: does Eigen's EigenSolver + expm_at_eigen agree
    // with LAPACK dgeev_ + the current column-scale+dgemm_ path, to
    // the plan's 3-decimal-place working tolerance? ---
    Matrix current_result = expm_at_current(lap.eigenvectors, lap.eigenvectors_inv, lap.eigenvalues_real, t);
    EMatd eigen_result_d = expm_at_eigen<double>(eig_vectors_d, eig_vectors_inv_d, eig_values_d, t);
    double diff_double = max_abs_diff(current_result, to_matrix(eigen_result_d));

    EMatf eigen_result_f = expm_at_eigen<float>(eig_vectors_f, eig_vectors_inv_f, eig_values_f, static_cast<float>(t));
    double diff_float = max_abs_diff(current_result, to_matrix(eigen_result_f.cast<double>()));

    std::cerr << "# correctness (max abs diff vs. current LAPACK-based result, t=" << t << "):\n";
    std::cerr << "#   Eigen double: " << diff_double << (diff_double < 5e-4 ? "  OK (<5e-4)" : "  FAIL") << "\n";
    std::cerr << "#   Eigen float:  " << diff_float << (diff_float < 5e-4 ? "  OK (<5e-4)" : "  FAIL") << "\n";

    // --- Timing ---
    std::cout
        << "operation,current_median_us,current_q1_us,current_q3_us,new_median_us,new_q1_us,new_q3_us,speedup\n";

    volatile double sink = 0;

    // (1) The dense multiply alone, current (dgemm_) vs. naive loop
    // vs. Eigen double vs. Eigen float.
    Stats mult_current = time_calls(
        [&]() {
            Matrix r = Matrix::mult(lap.eigenvectors, lap.eigenvectors_inv);
            sink += r(0, 0);
        },
        WARMUP, REPS);
    Stats mult_naive = time_calls(
        [&]() {
            Matrix r = naive_mult(lap.eigenvectors, lap.eigenvectors_inv);
            sink += r(0, 0);
        },
        WARMUP, REPS);
    Stats mult_eigen_d = time_calls(
        [&]() {
            EMatd r = eig_vectors_d * eig_vectors_inv_d;
            sink += r(0, 0);
        },
        WARMUP, REPS);
    Stats mult_eigen_f = time_calls(
        [&]() {
            EMatf r = eig_vectors_f * eig_vectors_inv_f;
            sink += r(0, 0);
        },
        WARMUP, REPS);
    report("dense_mult_20x20_naive_vs_current", mult_current, mult_naive);
    report("dense_mult_20x20_eigen_double_vs_current", mult_current, mult_eigen_d);
    report("dense_mult_20x20_eigen_float_vs_current", mult_current, mult_eigen_f);

    // (2) The full MatrixExpm::at(t)-equivalent operation: current vs.
    // Eigen double vs. Eigen float - this is what actually runs once
    // per Newton iteration per pair in likelihood_deriv().
    Stats expm_current = time_calls(
        [&]() {
            Matrix r = expm_at_current(lap.eigenvectors, lap.eigenvectors_inv, lap.eigenvalues_real, t);
            sink += r(0, 0);
        },
        WARMUP, REPS);
    Stats expm_eigen_d = time_calls(
        [&]() {
            EMatd r = expm_at_eigen<double>(eig_vectors_d, eig_vectors_inv_d, eig_values_d, t);
            sink += r(0, 0);
        },
        WARMUP, REPS);
    Stats expm_eigen_f = time_calls(
        [&]() {
            EMatf r = expm_at_eigen<float>(eig_vectors_f, eig_vectors_inv_f, eig_values_f, static_cast<float>(t));
            sink += r(0, 0);
        },
        WARMUP, REPS);
    report("expm_at_eigen_double_vs_current", expm_current, expm_eigen_d);
    report("expm_at_eigen_float_vs_current", expm_current, expm_eigen_f);

    if (sink == 123456.0)
    {
        std::cerr << "unreachable\n"; // prevent over-optimization
    }

    return 0;
}
