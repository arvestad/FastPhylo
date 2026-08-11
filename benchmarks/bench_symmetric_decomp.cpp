// Checks whether exploiting Q's time-reversibility (eq[i]*Q(i,j) ==
// eq[j]*Q(j,i), true for all 5 named models by construction; -F isn't
// implemented yet, see MatrixExpm's class comment) to decompose the
// symmetric similarity transform D^0.5 * Q * D^-0.5 (D=diag(eq)) via
// SelfAdjointEigenSolver - instead of the general, possibly-complex
// EigenSolver + an explicit .inverse() call - is actually cheaper in
// practice. Prompted by bench_decomp_share.cpp finding decomposition
// is 92-95% of total time at N=2, 20-26% at N=10, non-negligible even
// at N=25 - contrary to Q2's original "always negligible" finding,
// which was measured against the pre-optimization, much slower
// per-pair loop.
//
// This was the prototype (measured 3.15-3.8x). Implemented for real
// in MatrixExpm's second constructor (Matrix.hpp/.cpp, 2026-08-11) -
// kept here as the exploratory evidence, not the reference
// implementation.
//
// Standalone, exploratory - not wired into the CMake build. Build:
//   c++ -O2 -std=c++17 -I../include \
//     -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 \
//     bench_symmetric_decomp.cpp \
//     ../src/fastphylo/protein/Matrix.cpp \
//     ../src/fastphylo/protein/ModelMatrix.cpp \
//     -framework Accelerate -o bench_symmetric_decomp
//   ./bench_symmetric_decomp

#include <Eigen/Dense>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>
#include "fastphylo/protein/Matrix.hpp"
#include "fastphylo/protein/ModelMatrix.hpp"

using namespace std::chrono;

namespace
{

Eigen::Matrix<double, 20, 20> to_eigen(const Matrix &Q)
{
    Eigen::Matrix<double, 20, 20> m;
    for (int i = 0; i < 20; i++)
        for (int j = 0; j < 20; j++)
            m(i, j) = Q(i, j);
    return m;
}

double median_of(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

} // namespace

int main()
{
    const int WARMUP = 3;
    const int REPS = 51;

    struct
    {
        model_type type;
        const char *name;
    } models[] = {{wag, "wag"}, {jtt, "jtt"}, {day, "day"}, {mvr, "mvr"}, {lg, "lg"}};

    std::cout << "model,general_us,symmetric_us,speedup,max_abs_diff_expm_t1\n";

    for (auto &m : models)
    {
        Matrix Qraw = get_model_matrix(m.type);
        DblVec eq = get_model_vec(m.type);
        Eigen::Matrix<double, 20, 20> Q = to_eigen(Qraw);

        Eigen::Matrix<double, 20, 1> eqv;
        for (int i = 0; i < 20; i++)
            eqv(i) = eq[i];
        Eigen::Matrix<double, 20, 1> sqrt_eq = eqv.array().sqrt();
        Eigen::Matrix<double, 20, 20> D_half = sqrt_eq.asDiagonal();
        Eigen::Matrix<double, 20, 20> D_half_inv = sqrt_eq.array().inverse().matrix().asDiagonal();

        // --- current approach: general EigenSolver + explicit inverse ---
        for (int i = 0; i < WARMUP; i++)
        {
            Eigen::EigenSolver<Eigen::Matrix<double, 20, 20>> solver(Q);
            Eigen::Matrix<double, 20, 20> V = solver.eigenvectors().real();
            Eigen::Matrix<double, 20, 20> Vinv = V.inverse();
            (void)Vinv;
        }
        std::vector<double> general_samples;
        Eigen::Matrix<double, 20, 20> V_general, Vinv_general;
        Eigen::Matrix<double, 20, 1> eig_general;
        for (int i = 0; i < REPS; i++)
        {
            auto t0 = high_resolution_clock::now();
            Eigen::EigenSolver<Eigen::Matrix<double, 20, 20>> solver(Q);
            V_general = solver.eigenvectors().real();
            Vinv_general = V_general.inverse();
            eig_general = solver.eigenvalues().real();
            auto t1 = high_resolution_clock::now();
            general_samples.push_back(duration<double, std::micro>(t1 - t0).count());
        }

        // --- candidate: symmetric similarity transform + SelfAdjointEigenSolver ---
        for (int i = 0; i < WARMUP; i++)
        {
            Eigen::Matrix<double, 20, 20> S = D_half * Q * D_half_inv;
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 20, 20>> solver(S);
            (void)solver;
        }
        std::vector<double> sym_samples;
        Eigen::Matrix<double, 20, 20> V_sym, Vinv_sym;
        Eigen::Matrix<double, 20, 1> eig_sym;
        for (int i = 0; i < REPS; i++)
        {
            auto t0 = high_resolution_clock::now();
            Eigen::Matrix<double, 20, 20> S = D_half * Q * D_half_inv;
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 20, 20>> solver(S);
            // Q = D^-0.5 * S * D^0.5, S = U * L * U^T (U orthogonal)
            // => Q = (D^-0.5 * U) * L * (U^T * D^0.5) = V * L * V^-1
            V_sym = D_half_inv * solver.eigenvectors();
            Vinv_sym = solver.eigenvectors().transpose() * D_half;
            eig_sym = solver.eigenvalues();
            auto t1 = high_resolution_clock::now();
            sym_samples.push_back(duration<double, std::micro>(t1 - t0).count());
        }

        // Correctness check: compare exp(Q*1.0) from both decompositions.
        Eigen::Matrix<double, 20, 1> scale_g = (eig_general * 1.0).array().exp();
        Eigen::Matrix<double, 20, 20> expm_general = V_general * scale_g.asDiagonal() * Vinv_general;
        Eigen::Matrix<double, 20, 1> scale_s = (eig_sym * 1.0).array().exp();
        Eigen::Matrix<double, 20, 20> expm_sym = V_sym * scale_s.asDiagonal() * Vinv_sym;
        double max_abs_diff = (expm_general - expm_sym).cwiseAbs().maxCoeff();

        double gen_us = median_of(general_samples);
        double sym_us = median_of(sym_samples);
        std::cout << m.name << "," << gen_us << "," << sym_us << "," << (gen_us / sym_us) << "," << max_abs_diff
                  << "\n";
    }

    return 0;
}
