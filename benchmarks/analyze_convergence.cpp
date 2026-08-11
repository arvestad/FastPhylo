// Measures likelihood_calc()'s actual convergence behavior on real
// data: how many Newton iterations pairs typically need, how often
// MAX_ITERATIONS is hit, how often MIN_DISTANCE/MAX_DISTANCE clamp,
// and the achieved precision on t (|slope|/|curv| at the converged
// point, a first-order estimate of |t - t*| from the local quadratic
// model) - to check whether CONVERGENCE_TOL=0.1 is actually loose in
// practice, and how often MIN_DISTANCE=0.01 vs. a smaller floor would
// matter. Mirrors likelihood_calc()'s loop exactly (MaximumLikelihood.cpp)
// but calls the real likelihood_slope_curv() with iteration counting
// added, since likelihood_calc() itself doesn't expose that.
//
// Standalone, exploratory - not wired into the CMake build. Build:
//   c++ -O2 -std=c++17 -I../include -I<path-to-eigen3-include-dir> \
//     analyze_convergence.cpp ../src/fastphylo/protein/Matrix.cpp \
//     ../src/fastphylo/protein/ModelMatrix.cpp \
//     ../src/fastphylo/protein/MaximumLikelihood.cpp \
//     ../src/fastphylo/protein/ProtSeqCode.cpp \
//     ../src/fastphylo/protein/ProtSeqCompare.cpp \
//     -framework Accelerate -o analyze_convergence
//   ./analyze_convergence <fasta-file>

#include <algorithm>
#include <cfloat>
#include <cmath>
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

namespace
{

// Mirrors MaximumLikelihood.cpp's constants exactly (post the
// 2026-08-11 PAM-unit removal - substitutions/site throughout).
constexpr double MAX_DISTANCE = 5.0;
constexpr double MIN_DISTANCE = 0.001;
constexpr double CONVERGENCE_TOL = 0.1;
constexpr int MAX_ITERATIONS = 50;

double kimura_distance_local(double n_sum, double n_sum_diag)
{
    double d = (n_sum - n_sum_diag) / n_sum;
    double adjusted = d + (0.2 * pow(d, 2));
    adjusted = std::min(adjusted, 0.854);
    return -log(1 - adjusted);
}

enum class ExitReason
{
    converged,
    hit_min,
    hit_max,
    exhausted
};

struct Result
{
    ExitReason reason;
    int iterations;
    double t;
    double precision_estimate; // |slope|/|curv| at the returned point, for converged cases
};

// Exact mirror of likelihood_calc(), with iteration counting and
// precision-estimate tracking added.
Result likelihood_calc_instrumented(const Eigen::Matrix<float, 20, 20> &N, const MatrixExpm &Qdecomp)
{
    double n_sum = N.sum();
    double n_sum_diag = N.trace();
    if (n_sum - n_sum_diag < DBL_EPSILON)
    {
        return {ExitReason::converged, 0, 0, 0};
    }

    double t = kimura_distance_local(n_sum, n_sum_diag);
    if (t <= 0)
    {
        t = MIN_DISTANCE;
    }
    double delta = t / 2.0;

    for (int i = 0; i < MAX_ITERATIONS; i++)
    {
        LikelihoodDerivatives d = likelihood_slope_curv(N, Qdecomp, t);

        if (fabs(d.slope) < CONVERGENCE_TOL)
        {
            double precision = d.curv != 0 ? fabs(d.slope / d.curv) : 0;
            return {ExitReason::converged, i + 1, t, precision};
        }

        if (d.curv < 0)
        {
            t -= d.slope / d.curv;
            if (t > MAX_DISTANCE)
            {
                return {ExitReason::hit_max, i + 1, MAX_DISTANCE, 0};
            }
        }
        else
        {
            if ((d.slope > 0) != (delta > 0))
            {
                delta = -delta / 2;
            }
            t += delta;
        }

        if (t < MIN_DISTANCE)
        {
            return {ExitReason::hit_min, i + 1, MIN_DISTANCE, 0};
        }
    }
    return {ExitReason::exhausted, MAX_ITERATIONS, t, 0};
}

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

Matrix tally_to_matrix(const std::vector<std::size_t> &tally)
{
    Matrix m(20, 20);
    for (std::size_t a = 0; a < 20; a++)
        for (std::size_t b = 0; b < 20; b++)
            m(a, b) = static_cast<double>(tally[(a * 20) + b]);
    return m;
}

Eigen::Matrix<float, 20, 20> tally_to_eigen(const std::vector<std::size_t> &tally)
{
    Eigen::Matrix<float, 20, 20> m;
    for (std::size_t a = 0; a < 20; a++)
        for (std::size_t b = 0; b < 20; b++)
            m(a, b) = static_cast<float>(tally[(a * 20) + b]);
    return m;
}

template <typename T> void print_stats(std::vector<T> v, const char *name)
{
    if (v.empty())
    {
        std::cout << "  " << name << ": (none)\n";
        return;
    }
    std::sort(v.begin(), v.end());
    double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    std::cout << "  " << name << ": n=" << v.size() << " min=" << v.front() << " median=" << v[v.size() / 2]
              << " mean=" << mean << " p90=" << v[static_cast<std::size_t>(v.size() * 0.9)] << " max=" << v.back()
              << "\n";
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

    for (auto &m : models)
    {
        Matrix Q = get_model_matrix(m.type);
        DblVec eq = get_model_vec(m.type);
        MatrixExpm Qdecomp(Q, 1.0 / mean_substitution_rate(Q, eq));

        std::vector<int> iterations;
        std::vector<double> precisions;
        std::vector<double> min_distance_ts;
        int n_converged = 0, n_hit_min = 0, n_hit_max = 0, n_exhausted = 0, n_identical = 0;
        int pairs = 0;

        for (std::size_t i = 0; i < seqs.size(); i++)
        {
            for (std::size_t j = i + 1; j < seqs.size(); j++)
            {
                Eigen::Matrix<float, 20, 20> N = tally_to_eigen(ProtSeqCode::count_replacement_tally(
                    encoded[i].data(), encoded[i].size(), encoded[j].data(), encoded[j].size()));
                pairs++;
                Result r = likelihood_calc_instrumented(N, Qdecomp);
                if (r.iterations == 0)
                {
                    n_identical++;
                    continue;
                }
                iterations.push_back(r.iterations);
                switch (r.reason)
                {
                case ExitReason::converged:
                    n_converged++;
                    precisions.push_back(r.precision_estimate);
                    break;
                case ExitReason::hit_min:
                    n_hit_min++;
                    break;
                case ExitReason::hit_max:
                    n_hit_max++;
                    break;
                case ExitReason::exhausted:
                    n_exhausted++;
                    break;
                }
            }
        }

        std::cout << "=== " << m.name << " (" << pairs << " pairs, " << n_identical << " identical) ===\n";
        std::cout << "  converged=" << n_converged << " (" << (100.0 * n_converged / pairs) << "%)"
                  << "  hit_min=" << n_hit_min << " (" << (100.0 * n_hit_min / pairs) << "%)"
                  << "  hit_max=" << n_hit_max << " (" << (100.0 * n_hit_max / pairs) << "%)"
                  << "  exhausted=" << n_exhausted << " (" << (100.0 * n_exhausted / pairs) << "%)\n";
        print_stats(iterations, "iterations (all exit reasons)");
        print_stats(precisions, "precision estimate |slope/curv| at converged t (substitutions/site)");
    }

    return 0;
}
