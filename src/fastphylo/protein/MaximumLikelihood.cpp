#include "fastphylo/protein/MaximumLikelihood.hpp"
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iostream>
#include "fastphylo/protein/Matrix.hpp"

namespace
{
// "500 is infinity"/"1 is the smallest possible distance" are this
// codebase's pre-existing conventions (kimura_distance() below always
// predates this file); named here rather than left as magic numbers
// now that the loop that uses them is more than a few lines.
constexpr double MAX_DISTANCE = 500;
constexpr double MIN_DISTANCE = 1;
constexpr double CONVERGENCE_TOL = 0.001;
constexpr int MAX_ITERATIONS = 50;

// Mirrors dna_pairwise_sequence_likelihood.hpp's warnTooDiverged() -
// same UX principle (a saturated pair is a real, expected outcome on
// diverse real data, worth a stderr note rather than a silently
// clamped number), separate helper since this is a different module
// with its own distance scale/units.
void warn_protein_too_diverged()
{
    std::cerr << "warning: sequences too different for maximum-likelihood protein distance, using distance="
              << MAX_DISTANCE << std::endl;
}
} // namespace

/*!
 * Kimura distance, not corrected, used as a starting value for the
 * Newton search in likelihood_calc()
 * @param N The replacement count matrix, N(i,j) contains the number of actual
 *            replacements from amino acid i to amino acid j
 * @return The kimura distance in PAMs for the matrix N
 */
double kimura_distance(const Matrix &N)
{
    double n_sum = N.sum();
    double d = (n_sum - N.sum_diag()) / n_sum;

    double adjusted = d + (0.2 * pow(d, 2));
    adjusted = std::min(adjusted, 0.854); // Infinite distance

    adjusted = -100 * log(1 - adjusted);
    return adjusted;
}

/*!
 * Analytic first and second derivatives of the log-likelihood
 * log L(t) = sum_{i,j} N(i,j) log(P(t)(i,j)), where P(t) = e^{Qt}:
 *   P'(t)  = Q P(t)          (= P(t) Q - Q commutes with its own
 *                                matrix exponential)
 *   P''(t) = Q P'(t)
 *   slope  = d/dt   log L(t) = sum N .* P'(t)  ./ P(t)
 *   curv   = d^2/dt^2 log L(t) = sum N .* (P''(t) ./ P(t) - (P'(t) ./ P(t))^2)
 * Both derivatives are obtained from the same cached decomposition
 * (Qdecomp.at(t)) with no finite differencing - P''(t) costs one more
 * matrix multiply reusing P'(t), which slope alone already computed.
 * This replaces the previous version that only returned slope,
 * requiring likelihood_calc() to approximate curv by finite-
 * differencing slope at a fixed, tiny step - noise-dominated wherever
 * the true derivative is locally flat, which real diverged sequence
 * pairs commonly are. See planning/fastprot_ml_speedup_investigation_
 * plan.md's "Q3 results" for the real-data cases this was found on,
 * and PHYLIP's protdist.c (makedists()) for the reference approach
 * this follows.
 * @param N The replacement count matrix, N(i,j) contains the number of actual
 *            replacements from amino acid i to amino acid j
 * @param Q Rate matrix for replacements from one amino acid to another
 * @param Qdecomp Cached eigendecomposition of Q (see Matrix.hpp)
 * @param t The distance
 * @return The log-likelihood's slope and curvature at t
 */
LikelihoodDerivatives likelihood_slope_curv(const Matrix &N, const Matrix &Q, const MatrixExpm &Qdecomp, double t)
{
    Matrix pt = Qdecomp.at(t);
    Matrix p1t = Matrix::mult(pt, Q);
    Matrix p2t = Matrix::mult(p1t, Q);

    // Fused into one pass instead of elem_mult()/elem_div() (each
    // allocating a full temporary Matrix just to be summed and
    // discarded) - same per-element operation order as a direct
    // transcription of the formulas above.
    double slope = 0.0;
    double curv = 0.0;
    std::size_t n = N.get_rows() * N.get_cols();
    for (std::size_t k = 0; k < n; k++)
    {
        double n_k = N(k);
        if (n_k == 0.0)
        {
            continue; // contributes 0 to both sums either way - skip the division
        }
        // P(t) is a transition-probability matrix for a rate matrix
        // with strictly positive off-diagonal rates (every model this
        // is used with), so p > 0 for all t > 0 mathematically; floor
        // it defensively against floating-point underflow at extreme
        // t rather than letting a division produce inf/nan.
        double p = std::max(pt(k), DBL_MIN);
        double dp = p1t(k);
        double d2p = p2t(k);
        double ratio = dp / p;
        slope += n_k * ratio;
        curv += n_k * (d2p / p - (ratio * ratio));
    }
    return {slope, curv};
}

/*!
 * Computes the maximum-likelihood distance via a safeguarded
 * Newton-Raphson search against the analytic slope/curvature of the
 * log-likelihood (likelihood_slope_curv()), ported from PHYLIP's
 * protdist.c (makedists()).
 *
 * The previous version of this function approximated the log-
 * likelihood's second derivative by finite-differencing its first
 * derivative at a fixed, tiny step, and - if a single step's result
 * didn't look like an improvement by a shallow heuristic ("the
 * derivative is getting larger") - returned that one, unverified step
 * as final, without ever evaluating the function there. On real
 * (non-synthetic) data this was found to return distances wrong by
 * several tenths of a PAM-scaled unit for a large, model-dependent
 * fraction of pairs (see planning/fastprot_ml_speedup_investigation_
 * plan.md's "Q3 results").
 *
 * This version takes a true Newton step only where the curvature is
 * actually negative (i.e. genuinely concave, near a real maximum);
 * everywhere else (not yet close enough to the peak, or past an
 * inflection, where a quadratic model isn't trustworthy) it falls
 * back to a safe, direction-following step-halving search instead.
 * No step is ever trusted without evaluating the function there on
 * the next iteration - there is no "one unverified step and stop"
 * exit.
 * @param N The replacement count matrix, N(i,j) contains the number of actual
 *            replacements from amino acid i to amino acid j
 * @param Q Rate matrix for replacements from one amino acid to another
 * @return The distance with the maximum likelihood
 */
double likelihood_calc(const Matrix &N, const Matrix &Q, const MatrixExpm &Qdecomp)
{
    if (N.sum() - N.sum_diag() < DBL_EPSILON)
    {
        return 0;
    }

    double t = kimura_distance(N);
    if (t <= 0)
    {
        t = 1;
    }
    double delta = t / 2.0;

    for (int i = 0; i < MAX_ITERATIONS; i++)
    {
        LikelihoodDerivatives d = likelihood_slope_curv(N, Q, Qdecomp, t);

        if (fabs(d.slope) < CONVERGENCE_TOL)
        {
            return t;
        }

        if (d.curv < 0)
        {
            // Genuine local maximum nearby - the quadratic model is
            // trustworthy here, take the actual Newton step.
            t -= d.slope / d.curv;
            if (t > MAX_DISTANCE)
            {
                warn_protein_too_diverged();
                return MAX_DISTANCE;
            }
        }
        else
        {
            // Not concave here - a Newton step could go anywhere;
            // step in the direction slope indicates instead (log
            // likelihood increases with t iff slope > 0), halving
            // whenever the last step overshot (slope's sign now
            // disagrees with delta's direction).
            if ((d.slope > 0) != (delta > 0))
            {
                delta = -delta / 2;
            }
            t += delta;
        }

        if (t < MIN_DISTANCE)
        {
            return MIN_DISTANCE;
        }
    }
    return t;
}
