// How much does MatrixExpm(Q,eq)'s one-time decomposition cost as a
// share of calculate_ml_dists()'s total wall time, on the CURRENT
// (post all Eigen/float/redundancy optimizations) code - re-measured
// 2026-08-11 because Q2 (fastprot_ml_speedup_investigation_plan.md)
// measured this against the OLD, slower LAPACK-based per-pair loop;
// the per-pair loop has since gotten 2-3x faster, so decomposition's
// *share* of total time could have grown even though its absolute
// cost hasn't changed. Answers directly rather than trusting a stale
// percentage.
//
// Updated 2026-08-11 to use the symmetrized MatrixExpm(Q, eq, scale)
// constructor (bench_symmetric_decomp.cpp's prototype, now real -
// Matrix.hpp/.cpp) - this file's numbers now reflect current
// production decomposition cost, not the general-EigenSolver cost
// they originally measured. Separately, calculate_ml_dists() no
// longer builds its own decomposition at all (build_ml_decomposition()
// is now a distinct, explicitly-called step - ProtDistCalc.hpp), so a
// caller invoking it many times (bootstrap replicates, etc.) doesn't
// pay this per-call share repeatedly the way the numbers below might
// suggest - see fastprot_ml_speedup_implementation_plan.md's
// "build-once API" section for the real end-to-end measurement of
// that fix.
//
// Standalone, exploratory - not wired into the CMake build. Build:
//   c++ -O2 -std=c++17 -I../include \
//     -I/opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3 \
//     -I/opt/homebrew/Cellar/simde/0.8.2/include \
//     bench_decomp_share.cpp \
//     ../src/fastphylo/protein/Matrix.cpp \
//     ../src/fastphylo/protein/ModelMatrix.cpp \
//     ../src/fastphylo/protein/MaximumLikelihood.cpp \
//     ../src/fastphylo/protein/ProtSeqCode.cpp \
//     ../src/fastphylo/protein/ProtSeqCompare.cpp \
//     -framework Accelerate -o bench_decomp_share
//   ./bench_decomp_share <fasta-file>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include "fastphylo/protein/MaximumLikelihood.hpp"
#include "fastphylo/protein/Matrix.hpp"
#include "fastphylo/protein/ModelMatrix.hpp"
#include "fastphylo/protein/ProtSeqCode.hpp"
#include "fastphylo/protein/ProtSeqCompare.hpp"

using namespace std::chrono;

namespace
{

std::vector<std::string> read_fasta(const std::string &path)
{
    std::ifstream in(path);
    std::vector<std::string> seqs;
    std::string line, current;
    while (std::getline(in, line))
    {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.empty())
            continue;
        if (line[0] == '>')
        {
            if (!current.empty())
            {
                seqs.push_back(current);
                current.clear();
            }
        }
        else
            current += line;
    }
    if (!current.empty())
        seqs.push_back(current);
    return seqs;
}

Eigen::Matrix<float, 20, 20> tally_to_eigen(const std::vector<std::size_t> &tally)
{
    Eigen::Matrix<float, 20, 20> m;
    for (std::size_t a = 0; a < 20; a++)
        for (std::size_t b = 0; b < 20; b++)
            m(a, b) = static_cast<float>(tally[(a * 20) + b]);
    return m;
}

} // namespace

int main(int argc, char **argv)
{
    std::string path = argc >= 2 ? argv[1] : "../examples/globin_family.fasta";
    std::vector<std::string> seqs = read_fasta(path);
    std::cerr << "# loaded " << seqs.size() << " sequences from " << path << "\n";

    std::vector<std::vector<std::uint8_t>> encoded(seqs.size());
    for (std::size_t i = 0; i < seqs.size(); i++)
    {
        ProtSeqCode::encode_sequence(seqs[i], encoded[i]);
    }

    struct
    {
        model_type type;
        const char *name;
    } models[] = {{wag, "wag"}, {jtt, "jtt"}, {day, "day"}, {mvr, "mvr"}, {lg, "lg"}};

    const int DECOMP_REPS = 51;
    std::cout << "model,decomp_us,pairs,per_pair_us,total_pair_us,decomp_share_pct\n";

    for (auto &m : models)
    {
        Matrix Q = get_model_matrix(m.type);
        DblVec eq = get_model_vec(m.type);
        double scale = 1.0 / mean_substitution_rate(Q, eq);

        // Time decomposition alone, several reps.
        std::vector<double> decomp_samples;
        for (int i = 0; i < DECOMP_REPS; i++)
        {
            auto t0 = high_resolution_clock::now();
            MatrixExpm decomp(Q, eq, scale);
            auto t1 = high_resolution_clock::now();
            decomp_samples.push_back(duration<double, std::micro>(t1 - t0).count());
        }
        std::sort(decomp_samples.begin(), decomp_samples.end());
        double decomp_us = decomp_samples[decomp_samples.size() / 2];

        // Time the full per-pair loop (all Newton iterations, real data).
        MatrixExpm Qdecomp(Q, eq, scale);
        std::vector<Eigen::Matrix<float, 20, 20>> Ns;
        for (std::size_t i = 0; i < seqs.size(); i++)
            for (std::size_t j = i + 1; j < seqs.size(); j++)
                Ns.push_back(tally_to_eigen(ProtSeqCode::count_replacement_tally(
                    encoded[i].data(), encoded[i].size(), encoded[j].data(), encoded[j].size())));

        // Warmup
        double sink = 0;
        for (auto &N : Ns)
            sink += likelihood_calc(N, Qdecomp);

        auto t0 = high_resolution_clock::now();
        for (int rep = 0; rep < 5; rep++)
            for (auto &N : Ns)
                sink += likelihood_calc(N, Qdecomp);
        auto t1 = high_resolution_clock::now();
        double total_us = duration<double, std::micro>(t1 - t0).count() / 5.0;
        double per_pair_us = total_us / Ns.size();

        double decomp_share_pct = 100.0 * decomp_us / (decomp_us + total_us);
        std::cout << m.name << "," << decomp_us << "," << Ns.size() << "," << per_pair_us << "," << total_us << ","
                  << decomp_share_pct << "\n";
        if (sink == -1)
            std::cerr << ""; // prevent optimizing away
    }

    return 0;
}
