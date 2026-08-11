// Micro-experiment for Q3 (fastprot ML speedup investigation plan,
// planning/fastprot_ml_speedup_investigation_plan.md): the current
// finite-difference "Newton-Raphson" solver (likelihood_calc(),
// MaximumLikelihood.cpp) vs. Brent's method, both finding the root of
// likelihood_deriv(t) = 0 for real replacement-count matrices N.
//
// Standalone, exploratory - not wired into the CMake build. Build
// manually:
//   c++ -O2 -std=c++17 -I../include \
//     bench_ml_solver.cpp \
//     ../src/fastphylo/protein/Matrix.cpp \
//     ../src/fastphylo/protein/ModelMatrix.cpp \
//     ../src/fastphylo/protein/MaximumLikelihood.cpp \
//     ../src/fastphylo/protein/ProtSeqCode.cpp \
//     ../src/fastphylo/protein/ProtSeqCompare.cpp \
//     -framework Accelerate -o bench_ml_solver
//   ./bench_ml_solver <fasta-file> [model: wag|jtt|day|mvr|lg]

#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
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

// Minimal FASTA reader - sequences only, headers discarded.
// benchmarks/gen_dataset.py's output and examples/globin_family.fasta
// are both plain "one sequence per record" FASTA, so no multi-line
// wrapping handling is needed here.
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

// Mirrors ProtDistCalc.cpp's anonymous-namespace tally_to_matrix() -
// duplicated here (rather than exporting it from production code just
// for this experiment) since it's a few trivial lines.
Matrix tally_to_matrix(const std::vector<std::size_t> &tally)
{
    Matrix m(ProtSeqCode::NUM_CANONICAL_AA, ProtSeqCode::NUM_CANONICAL_AA);
    for (std::size_t a = 0; a < ProtSeqCode::NUM_CANONICAL_AA; a++)
    {
        for (std::size_t b = 0; b < ProtSeqCode::NUM_CANONICAL_AA; b++)
        {
            m(a, b) = static_cast<double>(tally[(a * ProtSeqCode::NUM_CANONICAL_AA) + b]);
        }
    }
    return m;
}

Matrix replacement_tally(const std::vector<std::uint8_t> &e1, const std::vector<std::uint8_t> &e2)
{
    return tally_to_matrix(ProtSeqCode::count_replacement_tally(e1.data(), e1.size(), e2.data(), e2.size()));
}

// A thin wrapper around likelihood_deriv() that counts calls - lets
// both solvers be compared on the metric that actually matters
// (evaluations of the expensive function), not just wall time, which
// is noisier at microsecond scale.
struct CountingDeriv
{
    const Matrix &N;
    const Matrix &Q;
    const MatrixExpm &Qdecomp;
    mutable long calls = 0;

    double operator()(double t) const
    {
        calls++;
        return likelihood_deriv(N, Q, Qdecomp, t);
    }
};

// Exactly mirrors likelihood_calc()'s Newton-Raphson loop
// (MaximumLikelihood.cpp), but calling through `f` so evaluations get
// counted, and returning whether it converged via the tol check vs.
// one of the three bail-out conditions.
enum class NewtonExit
{
    converged,
    hit_lower_bound,
    hit_upper_bound,
    deriv_growing,
    maxit_exhausted
};

struct NewtonResult
{
    double t;
    NewtonExit exit;
    long evals;
};

NewtonResult newton_solve(const CountingDeriv &f, double n_sum, double n_sum_diag)
{
    if (n_sum - n_sum_diag < DBL_EPSILON)
    {
        return {0, NewtonExit::converged, 0};
    }
    double d = (n_sum - n_sum_diag) / n_sum;
    double adjusted = d + (0.2 * pow(d, 2));
    adjusted = std::min(adjusted, 0.854);
    double t = -100 * log(1 - adjusted);
    if (t == 0)
    {
        t = 1;
    }

    double delta = 0.001;
    double tol = 0.001;
    int maxit = 50;

    double l_d = f(t);
    for (int i = 0; i < maxit; i++)
    {
        if (fabs(l_d) < tol)
        {
            return {t, NewtonExit::converged, f.calls};
        }
        double l_new = f(t + delta);
        double deriv = (l_new - l_d) / delta;
        t = t - (l_d / deriv);
        if (t < 1)
        {
            return {1, NewtonExit::hit_lower_bound, f.calls};
        }
        if (t > 500)
        {
            return {500, NewtonExit::hit_upper_bound, f.calls};
        }
        if (fabs(l_d) < fabs(l_new))
        {
            return {t, NewtonExit::deriv_growing, f.calls};
        }
        l_d = l_new;
    }
    return {t, NewtonExit::maxit_exhausted, f.calls};
}

// Classic Brent-Dekker root finder (see e.g. Numerical Recipes /
// Brent 1973), textbook form: combines bisection, secant, and inverse
// quadratic interpolation, with the bisection fallback guaranteeing
// convergence whenever [a,b] brackets a root. tol is on |b-a|,
// matching the Newton loop's use of tol as a convergence threshold
// (there it's on |f(t)|, here on the bracket width - the natural
// stopping criterion for a bracketing method).
struct BrentResult
{
    double t;
    bool bracket_valid;
    long evals;
    int iterations;
};

BrentResult brent_solve(const CountingDeriv &f, double a, double b, double tol, int maxit)
{
    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0)
    {
        return {b, false, f.calls, 0};
    }
    if (std::fabs(fa) < std::fabs(fb))
    {
        std::swap(a, b);
        std::swap(fa, fb);
    }
    double c = a;
    double fc = fa;
    bool mflag = true;
    double d = 0;
    int iter = 0;

    for (; iter < maxit; iter++)
    {
        if (fb == 0 || std::fabs(b - a) < tol)
        {
            break;
        }
        double s;
        if (fa != fc && fb != fc)
        {
            // Inverse quadratic interpolation
            s = (a * fb * fc) / ((fa - fb) * (fa - fc)) + (b * fa * fc) / ((fb - fa) * (fb - fc)) +
                (c * fa * fb) / ((fc - fa) * (fc - fb));
        }
        else
        {
            // Secant
            s = b - fb * (b - a) / (fb - fa);
        }
        double lo = (3 * a + b) / 4;
        double hi = b;
        if (lo > hi)
        {
            std::swap(lo, hi);
        }
        bool cond1 = !(s > lo && s < hi);
        bool cond2 = mflag && std::fabs(s - b) >= std::fabs(b - c) / 2;
        bool cond3 = !mflag && std::fabs(s - b) >= std::fabs(c - d) / 2;
        bool cond4 = mflag && std::fabs(b - c) < tol;
        bool cond5 = !mflag && std::fabs(c - d) < tol;
        if (cond1 || cond2 || cond3 || cond4 || cond5)
        {
            s = (a + b) / 2; // bisection
            mflag = true;
        }
        else
        {
            mflag = false;
        }
        double fs = f(s);
        d = c;
        c = b;
        fc = fb;
        if (fa * fs < 0)
        {
            b = s;
            fb = fs;
        }
        else
        {
            a = s;
            fa = fs;
        }
        if (std::fabs(fa) < std::fabs(fb))
        {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }
    return {b, true, f.calls, iter};
}

model_type parse_model(const std::string &s)
{
    if (s == "jtt")
        return jtt;
    if (s == "day")
        return day;
    if (s == "mvr")
        return mvr;
    if (s == "lg")
        return lg;
    return wag;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "usage: " << argv[0] << " <fasta-file> [model]\n";
        return 1;
    }
    std::string path = argv[1];
    model_type mt = argc >= 3 ? parse_model(argv[2]) : wag;

    std::vector<std::string> seqs = read_fasta(path);
    std::cerr << "# loaded " << seqs.size() << " sequences from " << path << "\n";

    std::vector<std::vector<std::uint8_t>> encoded(seqs.size());
    for (std::size_t i = 0; i < seqs.size(); i++)
    {
        ProtSeqCode::encode_sequence(seqs[i], encoded[i]);
    }

    Matrix Q = get_model_matrix(mt);
    MatrixExpm Qdecomp(Q);

    long total_newton_evals = 0;
    long total_brent_evals = 0;
    std::map<NewtonExit, int> newton_exit_counts;
    int bracket_failures = 0;
    int pairs = 0;
    double max_abs_diff = 0;
    double sum_wall_newton_us = 0;
    double sum_wall_brent_us = 0;
    std::size_t worst_i = 0;
    std::size_t worst_j = 0;
    Matrix worst_N(1, 1);

    for (std::size_t i = 0; i < seqs.size(); i++)
    {
        for (std::size_t j = i + 1; j < seqs.size(); j++)
        {
            Matrix N = replacement_tally(encoded[i], encoded[j]);
            pairs++;

            CountingDeriv fn{N, Q, Qdecomp};
            auto t0 = high_resolution_clock::now();
            NewtonResult nr = newton_solve(fn, N.sum(), N.sum_diag());
            auto t1 = high_resolution_clock::now();
            sum_wall_newton_us += duration<double, std::micro>(t1 - t0).count();
            total_newton_evals += nr.evals;
            newton_exit_counts[nr.exit]++;

            CountingDeriv fb{N, Q, Qdecomp};
            double tol = 0.001;
            auto t2 = high_resolution_clock::now();
            BrentResult br = brent_solve(fb, 1.0, 500.0, tol, 100);
            auto t3 = high_resolution_clock::now();
            sum_wall_brent_us += duration<double, std::micro>(t3 - t2).count();
            total_brent_evals += br.evals;
            if (!br.bracket_valid)
            {
                bracket_failures++;
                continue; // no meaningful t to compare
            }

            // Compare at the *reported* distance scale (0.01x, per
            // calculate_ml_dists()'s scaling) to the 3-decimal-place
            // working tolerance.
            double newton_dist = 0.01 * nr.t;
            double brent_dist = 0.01 * br.t;
            double diff = std::fabs(newton_dist - brent_dist);
            if (diff >= max_abs_diff)
            {
                max_abs_diff = diff;
                worst_i = i;
                worst_j = j;
                worst_N = N;
            }
        }
    }

    std::cout << "pairs," << pairs << "\n";
    std::cout << "bracket_failures," << bracket_failures << " (" << (100.0 * bracket_failures / pairs) << "%)\n";
    std::cout << "newton_total_evals," << total_newton_evals << "\n";
    std::cout << "brent_total_evals," << total_brent_evals << "\n";
    std::cout << "newton_avg_evals_per_pair," << (double)total_newton_evals / pairs << "\n";
    std::cout << "brent_avg_evals_per_pair," << (double)total_brent_evals / pairs << "\n";
    std::cout << "newton_total_wall_us," << sum_wall_newton_us << "\n";
    std::cout << "brent_total_wall_us," << sum_wall_brent_us << "\n";
    std::cout << "max_abs_distance_diff," << max_abs_diff
              << (max_abs_diff < 5e-4 ? "  OK (<5e-4, 3-decimal tolerance)" : "  EXCEEDS TOLERANCE") << "\n";
    std::cout << "newton_exit_reasons:\n";
    const char *names[] = {"converged", "hit_lower_bound(t<1)", "hit_upper_bound(t>500)", "deriv_growing(bailout)",
                           "maxit_exhausted"};
    for (auto &[exit, count] : newton_exit_counts)
    {
        std::cout << "  " << names[static_cast<int>(exit)] << "," << count << " (" << (100.0 * count / pairs) << "%)\n";
    }

    // Diagnostic: re-solve and dump the likelihood_deriv(t) curve for
    // the single worst-disagreement pair, to see whether Newton's
    // bail-out landed somewhere genuinely different from the function's
    // actual root, or whether the two solvers are both "reasonable"
    // given a badly-conditioned or non-monotonic derivative.
    if (worst_N.get_rows() > 1)
    {
        std::cout << "\n# worst-disagreement pair: seq[" << worst_i << "] vs seq[" << worst_j
                  << "], diff=" << max_abs_diff << "\n";
        std::cout << "# N.sum()=" << worst_N.sum() << " N.sum_diag()=" << worst_N.sum_diag()
                  << " (observed substitution fraction=" << (1.0 - worst_N.sum_diag() / worst_N.sum()) << ")\n";
        CountingDeriv fn{worst_N, Q, Qdecomp};
        NewtonResult nr = newton_solve(fn, worst_N.sum(), worst_N.sum_diag());
        CountingDeriv fb{worst_N, Q, Qdecomp};
        BrentResult br = brent_solve(fb, 1.0, 500.0, 0.001, 100);
        std::cout << "# newton: t=" << nr.t << " exit=" << names[static_cast<int>(nr.exit)] << " evals=" << nr.evals
                  << "\n";
        std::cout << "# brent:  t=" << br.t << " evals=" << br.evals << " iterations=" << br.iterations << "\n";
        std::cout << "# likelihood_deriv(t) grid:\n";
        std::cout << "t,likelihood_deriv\n";
        for (double t : {1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 150.0, 200.0, 300.0, 400.0, 500.0})
        {
            std::cout << t << "," << likelihood_deriv(worst_N, Q, Qdecomp, t) << "\n";
        }
    }

    return 0;
}
