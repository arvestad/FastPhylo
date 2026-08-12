// Regression test for a real, reproducible bug (2026-08-12): fastdist's
// default K2P distance (fixed transition/transversion ratio, via -T,
// default 2.0) used an unbounded secant search that could overshoot far
// enough for exp(2*x*K) to overflow, either producing nan or - worse -
// a spurious exact-zero derivative from inf-vs-inf cancellation that the
// convergence check mistook for a true root, silently returning a
// wildly wrong finite distance (e.g. ~56 instead of ~0.1). Both failure
// modes were 100% reproducible via the real `fastdist -T <ratio far
// from the data>` CLI, not just synthetic unit inputs - see
// Kimura2parameter.cpp's secant_search() for the fix (bounded search,
// falls back to TOO_DIVERGED_DISTANCE/warnTooDiverged() like every
// other distance formula in that file already does for its own
// degenerate cases).
#include "fastphylo/dna/dna_pairwise_sequence_likelihood.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>

static void test_default_ratio_sane()
{
    // tstvratio=2 (K=4, fastdist's default) on a realistic, moderately
    // diverged pair - must stay far below TOO_DIVERGED_DISTANCE and be
    // a normal, finite, non-negative value. This is the fix's "don't
    // regress the common case" checkpoint.
    simple_string_distance sd{0, 154.0F, 27.0F}; // ts=154, tv=27, n implied by strlen below
    ML_string_distance tp = compute_K2P_fixratio(2000, sd, 2.0F);
    assert(std::isfinite(tp.distance));
    assert(tp.distance > 0.0 && tp.distance < 1.0);
}

static void test_mismatched_ratio_reports_too_diverged_not_garbage()
{
    // The exact real-world failure case: a fixed ratio (tstvratio=8,
    // K=16) that doesn't match the data's actual transition/
    // transversion split (ts=21, tv=161 - transversion-heavy, the
    // opposite of what K=16 assumes). Before the fix: nan. Must now be
    // a finite, bounded, honest TOO_DIVERGED_DISTANCE - not garbage,
    // and not silently wrong.
    simple_string_distance sd{0, 21.0F, 161.0F};
    ML_string_distance tp =
        compute_K2P_fixratio(2000, sd, 8.0F); // fixRatio=8 -> K=16, matches the `fastdist -T 8` reproduction
    assert(std::isfinite(tp.distance));
    assert(tp.distance <= TOO_DIVERGED_DISTANCE + 1e-6);
}

static void test_false_convergence_case_now_bounded()
{
    // The second, more dangerous failure mode found during triage: a
    // secant overshoot that landed exactly where floating-point
    // overflow made the derivative evaluate to 0, which the old,
    // unbounded loop mistook for a converged true root - silently
    // returning distance ~55.87 for a pair whose real K2P distance is
    // ~0.18. fixRatio=0.125 -> K=0.25; ts=204, tv=115.
    simple_string_distance sd{0, 204.0F, 115.0F};
    ML_string_distance tp = compute_K2P_fixratio(2000, sd, 0.125F);
    assert(std::isfinite(tp.distance));
    assert(tp.distance <= TOO_DIVERGED_DISTANCE + 1e-6);
}

static void test_replacement_probabilities_stay_valid_when_too_diverged()
{
    // Even when the distance itself saturates to TOO_DIVERGED_DISTANCE,
    // the change-observation-probability fields (used elsewhere, e.g.
    // ambiguity resolution) must stay well-formed probabilities, not
    // nan/inf/out-of-range - Kimura2parameter.cpp's own asserts already
    // guard this at runtime; this test just exercises that path
    // explicitly under ctest rather than relying on incidentally
    // hitting it.
    simple_string_distance sd{0, 21.0F, 161.0F};
    ML_string_distance tp = compute_K2P_fixratio(2000, sd, 8.0F);
    auto valid_prob = [](float p) { return std::isfinite(p) && p >= 0.0F && p <= 1.0F; };
    assert(valid_prob(tp.A_A));
    assert(valid_prob(tp.A_C));
    assert(valid_prob(tp.A_G));
    assert(valid_prob(tp.A_T));
}

int main()
{
    test_default_ratio_sane();
    test_mismatched_ratio_reports_too_diverged_not_garbage();
    test_false_convergence_case_now_bounded();
    test_replacement_probabilities_stay_valid_when_too_diverged();
    printf("Kimura2parameter_test: all tests passed\n");
    return 0;
}
