// Re-baseline microbenchmark for the Eigen matrix-backend
// implementation (planning/fastprot_ml_speedup_implementation_plan.md):
// times the CURRENT, post-solver-fix likelihood_slope_curv()
// (MaximumLikelihood.cpp - now three matrix operations: P(t), P'(t),
// P''(t), not the two Q1's original microbenchmark measured) against
// an Eigen-based reimplementation of the same function, both double
// and float, on real data. Confirms (or corrects) the extrapolated
// 1.2-1.4x/2.4-2.8x estimate before committing to the float step.
//
// Standalone, exploratory - not wired into the CMake build. Build
// manually:
//   c++ -O2 -std=c++17 -I../include \
//     -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 \
//     bench_ml_hotpath.cpp \
//     ../src/fastphylo/protein/Matrix.cpp \
//     ../src/fastphylo/protein/ModelMatrix.cpp \
//     ../src/fastphylo/protein/MaximumLikelihood.cpp \
//     ../src/fastphylo/protein/ProtSeqCode.cpp \
//     ../src/fastphylo/protein/ProtSeqCompare.cpp \
//     -framework Accelerate -o bench_ml_hotpath
//   ./bench_ml_hotpath

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "fastphylo/protein/MaximumLikelihood.hpp"
#include "fastphylo/protein/Matrix.hpp"
#include "fastphylo/protein/ModelMatrix.hpp"
#include "fastphylo/protein/ProtSeqCode.hpp"
#include "fastphylo/protein/ProtSeqCompare.hpp"

using namespace std::chrono;

namespace
{

constexpr std::size_t N_AA = 20;
template <typename Scalar> using EMat = Eigen::Matrix<Scalar, N_AA, N_AA>;
template <typename Scalar> using EVec = Eigen::Matrix<Scalar, N_AA, 1>;

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

EMat<double> to_eigen(const Matrix &m)
{
    EMat<double> e;
    for (std::size_t r = 0; r < N_AA; r++)
    {
        for (std::size_t c = 0; c < N_AA; c++)
        {
            e(static_cast<int>(r), static_cast<int>(c)) = m(r, c);
        }
    }
    return e;
}

// Eigen-based equivalent of likelihood_slope_curv() - mirrors
// MaximumLikelihood.cpp's production function exactly (same
// multiplication order, same fused summation), but with
// Eigen::Matrix<Scalar,20,20> internally instead of Matrix/dgemm_.
template <typename Scalar>
LikelihoodDerivatives eigen_slope_curv(const Matrix &N, const EMat<Scalar> &eigenvectors,
                                        const EMat<Scalar> &eigenvectors_inv, const EVec<Scalar> &eigenvalues,
                                        const EMat<Scalar> &Q_eigen, Scalar t)
{
    EVec<Scalar> scale = (eigenvalues * t).array().exp();
    EMat<Scalar> pt = eigenvectors * scale.asDiagonal() * eigenvectors_inv;
    EMat<Scalar> p1t = pt * Q_eigen;
    EMat<Scalar> p2t = p1t * Q_eigen;

    double slope = 0.0;
    double curv = 0.0;
    for (int r = 0; r < static_cast<int>(N_AA); r++)
    {
        for (int c = 0; c < static_cast<int>(N_AA); c++)
        {
            double n_k = N(r, c);
            if (n_k == 0.0)
            {
                continue;
            }
            double p = std::max(static_cast<double>(pt(r, c)), 1e-300);
            double dp = static_cast<double>(p1t(r, c));
            double d2p = static_cast<double>(p2t(r, c));
            double ratio = dp / p;
            slope += n_k * ratio;
            curv += n_k * (d2p / p - (ratio * ratio));
        }
    }
    return {slope, curv};
}

std::vector<std::string> read_fasta(const std::string &path)
{
    std::ifstream in(path);
    std::vector<std::string> seqs;
    std::string line;
    std::string current;
    while (std::getline(in, line))
    {
        if (!line.empty() && line.back() == '\r')
        {
            line.pop_back();
        }
        if (line.empty())
        {
            continue;
        }
        if (line[0] == '>')
        {
            if (!current.empty())
            {
                seqs.push_back(current);
                current.clear();
            }
        }
        else
        {
            current += line;
        }
    }
    if (!current.empty())
    {
        seqs.push_back(current);
    }
    return seqs;
}

Matrix tally_to_matrix(const std::vector<std::size_t> &tally)
{
    Matrix m(N_AA, N_AA);
    for (std::size_t a = 0; a < N_AA; a++)
    {
        for (std::size_t b = 0; b < N_AA; b++)
        {
            m(a, b) = static_cast<double>(tally[(a * N_AA) + b]);
        }
    }
    return m;
}

} // namespace

int main(int argc, char **argv)
{
    std::string path = argc >= 2 ? argv[1] : "../examples/globin_family.fasta";
    std::vector<std::string> seqs = read_fasta(path);
    if (seqs.size() < 2)
    {
        std::cerr << "need at least 2 sequences in " << path << "\n";
        return 1;
    }
    std::vector<std::uint8_t> e0;
    std::vector<std::uint8_t> e1;
    ProtSeqCode::encode_sequence(seqs[0], e0);
    ProtSeqCode::encode_sequence(seqs[1], e1);
    Matrix N = tally_to_matrix(ProtSeqCode::count_replacement_tally(e0.data(), e0.size(), e1.data(), e1.size()));

    const int WARMUP = 3;
    const int REPS = 21;
    const double t = 100.0; // representative mid-range branch length

    struct
    {
        model_type type;
        const char *name;
    } models[] = {{wag, "wag"}, {jtt, "jtt"}, {day, "day"}, {mvr, "mvr"}, {lg, "lg"}};

    std::cout << "model,current_median_us,current_q1_us,current_q3_us,eigen_median_us,eigen_q1_us,eigen_q3_us,"
                 "speedup\n";

    for (auto &m : models)
    {
        Matrix Q = get_model_matrix(m.type);
        MatrixExpm Qdecomp(Q);
        EMat<double> Q_eigen_d = to_eigen(Q);

        Eigen::EigenSolver<EMat<double>> es(Q_eigen_d);
        EMat<double> eigvec_d = es.eigenvectors().real();
        EMat<double> eigvec_inv_d = eigvec_d.inverse();
        EVec<double> eigval_d = es.eigenvalues().real();

        EMat<float> Q_eigen_f = Q_eigen_d.cast<float>();
        EMat<float> eigvec_f = eigvec_d.cast<float>();
        EMat<float> eigvec_inv_f = eigvec_inv_d.cast<float>();
        EVec<float> eigval_f = eigval_d.cast<float>();

        // Correctness check against the production (LAPACK-based) result.
        LikelihoodDerivatives current = likelihood_slope_curv(N, Q, Qdecomp, t);
        LikelihoodDerivatives eig_d = eigen_slope_curv<double>(N, eigvec_d, eigvec_inv_d, eigval_d, Q_eigen_d, t);
        LikelihoodDerivatives eig_f = eigen_slope_curv<float>(N, eigvec_f, eigvec_inv_f, eigval_f, Q_eigen_f,
                                                               static_cast<float>(t));
        double diff_d = std::max(std::fabs(current.slope - eig_d.slope), std::fabs(current.curv - eig_d.curv));
        double diff_f = std::max(std::fabs(current.slope - eig_f.slope), std::fabs(current.curv - eig_f.curv));

        volatile double sink = 0;
        Stats cur_stats = time_calls(
            [&]() {
                LikelihoodDerivatives d = likelihood_slope_curv(N, Q, Qdecomp, t);
                sink += d.slope;
            },
            WARMUP, REPS);
        Stats eig_d_stats = time_calls(
            [&]() {
                LikelihoodDerivatives d = eigen_slope_curv<double>(N, eigvec_d, eigvec_inv_d, eigval_d, Q_eigen_d, t);
                sink += d.slope;
            },
            WARMUP, REPS);
        Stats eig_f_stats = time_calls(
            [&]() {
                LikelihoodDerivatives d = eigen_slope_curv<float>(N, eigvec_f, eigvec_inv_f, eigval_f, Q_eigen_f,
                                                                   static_cast<float>(t));
                sink += d.slope;
            },
            WARMUP, REPS);

        report(std::string(m.name) + "_double", cur_stats, eig_d_stats);
        std::cout << "  # max|slope,curv diff| double=" << diff_d << "\n";
        report(std::string(m.name) + "_float", cur_stats, eig_f_stats);
        std::cout << "  # max|slope,curv diff| float=" << diff_f << "\n";

        if (sink == 123456.0)
        {
            std::cerr << "unreachable\n";
        }
    }

    return 0;
}
