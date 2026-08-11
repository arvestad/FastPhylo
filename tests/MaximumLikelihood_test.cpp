// Regression test for likelihood_calc()'s safeguarded Newton-Raphson
// search (MaximumLikelihood.cpp).
//
// The property under test: likelihood_calc() must never return a
// distance it hasn't actually verified - either the log-likelihood's
// slope is genuinely near zero there (a real root), or the answer is
// a legitimate boundary clamp (MIN_DISTANCE/MAX_DISTANCE). The
// previous implementation violated this on real data: a "derivative
// is getting larger" heuristic could return a single, unevaluated
// Newton step as final. This test runs likelihood_calc() on every
// pair of a real protein family (not the narrow-divergence synthetic
// benchmarks that never triggered the bug) across all 5 ML models,
// and independently checks the returned distance against that
// invariant - see planning/fastprot_ml_speedup_investigation_plan.md's
// "Q3 results" for the investigation this came out of.
//
// GLOBIN_FIXTURE_PATH is supplied by CMake (see src/c++/CMakeLists.txt)
// as an absolute path to examples/globin_family.fasta, so this test
// doesn't depend on ctest's working directory.

#undef NDEBUG
#include <cassert>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "fastphylo/protein/MaximumLikelihood.hpp"
#include "fastphylo/protein/Matrix.hpp"
#include "fastphylo/protein/ModelMatrix.hpp"
#include "fastphylo/protein/ProtDistCalc.hpp"
#include "fastphylo/protein/ProtSeqCode.hpp"
#include "fastphylo/protein/ProtSeqCompare.hpp"

namespace
{

// Mirrors MaximumLikelihood.cpp's own constants (substitutions/site).
constexpr double MAX_DISTANCE = 5.0;
constexpr double MIN_DISTANCE = 0.001;
constexpr double CONVERGENCE_TOL = 0.1;
// A returned t is accepted as "at" a boundary if it's within this of
// it - likelihood_calc() returns MIN_DISTANCE/MAX_DISTANCE exactly
// today, but this test's job is to check the invariant, not pin down
// exact floating-point equality to an implementation detail.
constexpr double BOUNDARY_SLACK = 1e-9;

// Minimal FASTA reader - sequences only, headers discarded. Mirrors
// benchmarks/gen_dataset.py's/examples/globin_family.fasta's plain
// "one sequence per record" format.
std::vector<std::string> read_fasta(const std::string &path)
{
    std::ifstream in(path);
    assert(in.is_open() && "fixture file must exist - see GLOBIN_FIXTURE_PATH");
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

// Mirrors ProtDistCalc.cpp's anonymous-namespace tally_to_eigen() -
// duplicated here rather than exporting it from production code just
// for this test, same reasoning as ProtSeqCompare_test.cpp's oracle.
Eigen::Matrix<float, 20, 20> tally_to_eigen(const std::vector<std::size_t> &tally)
{
    Eigen::Matrix<float, 20, 20> m;
    for (std::size_t a = 0; a < ProtSeqCode::NUM_CANONICAL_AA; a++)
    {
        for (std::size_t b = 0; b < ProtSeqCode::NUM_CANONICAL_AA; b++)
        {
            m(a, b) = static_cast<float>(tally[(a * ProtSeqCode::NUM_CANONICAL_AA) + b]);
        }
    }
    return m;
}

Eigen::Matrix<float, 20, 20> replacement_tally(const std::vector<std::uint8_t> &e1, const std::vector<std::uint8_t> &e2)
{
    return tally_to_eigen(ProtSeqCode::count_replacement_tally(e1.data(), e1.size(), e2.data(), e2.size()));
}

// The property under test, isolated so both the "does it hold" check
// and its failure message are in one place.
void assert_verified_answer(const Eigen::Matrix<float, 20, 20> &N, const MatrixExpm &Qdecomp, double t,
                            const char *model_name, std::size_t i, std::size_t j)
{
    if (N.sum() - N.trace() < DBL_EPSILON)
    {
        assert(t == 0 && "identical sequences must return distance 0");
        return;
    }
    bool at_min = std::fabs(t - MIN_DISTANCE) < BOUNDARY_SLACK;
    bool at_max = std::fabs(t - MAX_DISTANCE) < BOUNDARY_SLACK;
    if (at_min || at_max)
    {
        return; // a legitimate, disclosed boundary clamp - not a silent guess
    }
    LikelihoodDerivatives d = likelihood_slope_curv(N, Qdecomp, t);
    if (std::fabs(d.slope) >= CONVERGENCE_TOL)
    {
        std::cerr << "FAIL: model=" << model_name << " pair=(" << i << "," << j << ") t=" << t << " slope=" << d.slope
                  << " (expected |slope| < " << CONVERGENCE_TOL
                  << ", or t at a boundary) - likelihood_calc() returned an unverified answer\n";
        assert(false);
    }
}

void test_globin_family_all_models()
{
    std::string path = GLOBIN_FIXTURE_PATH;
    std::vector<std::string> seqs = read_fasta(path);
    assert(seqs.size() > 1);

    std::vector<std::vector<std::uint8_t>> encoded(seqs.size());
    for (std::size_t i = 0; i < seqs.size(); i++)
    {
        ProtSeqCode::encode_sequence(seqs[i], encoded[i]);
    }

    struct
    {
        model_type type;
        const char *name;
    } models[] = {{wag, "wag"},
                  {jtt, "jtt"},
                  {day, "day"},
                  {mvr, "mvr"},
                  {lg, "lg"},
                  {jtt_dcmut, "jtt_dcmut"},
                  {vt, "vt"},
                  {hivb, "hivb"},
                  {hivw, "hivw"},
                  {cprev, "cprev"},
                  {blosum62, "blosum62"},
                  {dcmut, "dcmut"},
                  {mtrev, "mtrev"},
                  {rtrev, "rtrev"},
                  {pmb, "pmb"}};

    for (auto &m : models)
    {
        MatrixExpm Qdecomp = build_ml_decomposition(m.type);
        std::size_t checked = 0;
        for (std::size_t i = 0; i < seqs.size(); i++)
        {
            for (std::size_t j = i + 1; j < seqs.size(); j++)
            {
                Eigen::Matrix<float, 20, 20> N = replacement_tally(encoded[i], encoded[j]);
                double t = likelihood_calc(N, Qdecomp);
                assert_verified_answer(N, Qdecomp, t, m.name, i, j);
                checked++;
            }
        }
        std::cerr << "  " << m.name << ": " << checked << " pairs verified\n";
    }
}

// Direct sanity check on kimura_distance()/likelihood_calc()'s
// identical-sequence short-circuit, independent of fixture data.
void test_identical_sequences_give_zero_distance()
{
    MatrixExpm Qdecomp = build_ml_decomposition(wag);
    Eigen::Matrix<float, 20, 20> N = Eigen::Matrix<float, 20, 20>::Zero();
    for (std::size_t a = 0; a < ProtSeqCode::NUM_CANONICAL_AA; a++)
    {
        N(a, a) = 42; // every position matches, no replacements observed
    }
    assert(likelihood_calc(N, Qdecomp) == 0);
}

} // namespace

int main()
{
    test_identical_sequences_give_zero_distance();
    test_globin_family_all_models();
    std::cerr << "MaximumLikelihood_test: all checks passed\n";
    return 0;
}
