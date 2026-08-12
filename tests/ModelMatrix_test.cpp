// Regression test for the embedded protein model matrices
// (ModelMatrix.cpp): every named model's rate matrix Q must have each
// row sum to zero (Q(i,i) == -sum_{j!=i} Q(i,j)) - the fundamental
// invariant of a valid continuous-time Markov rate matrix - and its
// equilibrium frequency vector must sum to ~1.
//
// Added after finding JTT's Q(0,0) was hardcoded as -1.24783051471057305,
// exactly 100x the correct value (-0.01247830514710573, the actual sum
// of row 0's off-diagonal entries) - a decimal-point transcription
// error that silently corrupted every JTT maximum-likelihood distance
// (cross-checked and fixed against both PHYLIP's protdist and
// IQ-TREE2's model/modelprotein.cpp). No existing test caught this:
// MaximumLikelihood_test.cpp only checks the solver's self-consistency
// (does the returned answer match the model's own slope=0), which
// holds even for a wrong model. This test checks the model data
// itself, independent of the solver.

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <iostream>
#include <numeric>
#include "fastphylo/protein/Matrix.hpp"
#include "fastphylo/protein/ModelMatrix.hpp"
#include "fastphylo/protein/ProtSeqCode.hpp"

namespace
{

void test_rows_sum_to_zero()
{
    struct
    {
        model_type type;
        const char *name;
    } models[] = {{wag, "WAG"},
                  {jtt, "JTT"},
                  {day, "Dayhoff"},
                  {mvr, "MVR"},
                  {lg, "LG"},
                  {jtt_dcmut, "JTT-DCMUT"},
                  {vt, "VT"},
                  {hivb, "HIVb"},
                  {hivw, "HIVw"},
                  {cprev, "cpREV"},
                  {blosum62, "BLOSUM62"},
                  {dcmut, "DCMUT"},
                  {mtrev, "MtREV"},
                  {rtrev, "RtREV"},
                  {pmb, "PMB"}};

    for (auto &m : models)
    {
        Matrix Q = get_model_matrix(m.type);
        for (std::size_t row = 0; row < ProtSeqCode::NUM_CANONICAL_AA; row++)
        {
            double off_diag_sum = 0;
            for (std::size_t col = 0; col < ProtSeqCode::NUM_CANONICAL_AA; col++)
            {
                if (col != row)
                {
                    off_diag_sum += Q(row, col);
                }
            }
            double residual = Q(row, row) + off_diag_sum;
            if (std::fabs(residual) >= 1e-6)
            {
                std::cerr << "FAIL: model=" << m.name << " row=" << row << " diag=" << Q(row, row)
                          << " off_diag_sum=" << off_diag_sum << " residual=" << residual
                          << " (expected ~0 - row does not sum to zero)\n";
                assert(false);
            }
        }
        std::cerr << "  " << m.name << ": all rows sum to zero\n";
    }
}

void test_equilibrium_frequencies_sum_to_one()
{
    struct
    {
        model_type type;
        const char *name;
    } models[] = {{wag, "WAG"},
                  {jtt, "JTT"},
                  {day, "Dayhoff"},
                  {mvr, "MVR"},
                  {lg, "LG"},
                  {jtt_dcmut, "JTT-DCMUT"},
                  {vt, "VT"},
                  {hivb, "HIVb"},
                  {hivw, "HIVw"},
                  {cprev, "cpREV"},
                  {blosum62, "BLOSUM62"},
                  {dcmut, "DCMUT"},
                  {mtrev, "MtREV"},
                  {rtrev, "RtREV"},
                  {pmb, "PMB"}};

    for (auto &m : models)
    {
        DblVec eq = get_model_vec(m.type);
        assert(eq.size() == ProtSeqCode::NUM_CANONICAL_AA);
        double sum = std::accumulate(eq.begin(), eq.end(), 0.0);
        if (std::fabs(sum - 1.0) >= 1e-2)
        {
            std::cerr << "FAIL: model=" << m.name << " equilibrium frequencies sum to " << sum << " (expected ~1)\n";
            assert(false);
        }
    }
    std::cerr << "  all models: equilibrium frequencies sum to ~1\n";
}

} // namespace

int main()
{
    test_rows_sum_to_zero();
    test_equilibrium_frequencies_sum_to_one();
    std::cerr << "ModelMatrix_test: all checks passed\n";
    return 0;
}
