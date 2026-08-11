#include "fastphylo/protein/MaximumLikelihood.hpp"
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <limits>
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
 * Newton search in likelihood_calc(). Takes n_sum/n_sum_diag rather
 * than the replacement-count matrix itself - likelihood_calc() always
 * computes both first for its own early-exit check, so recomputing
 * them here from N again (a former version of this function did) was
 * pure duplicated work, once per pair.
 * @param n_sum Sum of all entries of the replacement count matrix N
 * @param n_sum_diag Sum of N's diagonal entries
 * @return The kimura distance in PAMs
 */
double kimura_distance(double n_sum, double n_sum_diag)
{
    double d = (n_sum - n_sum_diag) / n_sum;

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
 *   P''(t) = Q P'(t) = P(t) Q^2 (computed this way - directly from
 *                                P(t), not chained through P'(t) - so
 *                                P'(t)/P''(t) don't depend on each
 *                                other and can be evaluated
 *                                independently; Q^2 is precomputed
 *                                once per run, see Matrix.hpp)
 *   slope  = d/dt   log L(t) = sum N .* P'(t)  ./ P(t)
 *   curv   = d^2/dt^2 log L(t) = sum N .* (P''(t) ./ P(t) - (P'(t) ./ P(t))^2)
 * Both derivatives are obtained from the same cached decomposition
 * with no finite differencing. This replaces a version that only
 * returned slope, requiring likelihood_calc() to approximate curv by
 * finite-differencing slope at a fixed, tiny step - noise-dominated
 * wherever the true derivative is locally flat, which real diverged
 * sequence pairs commonly are. See planning/fastprot_ml_speedup_
 * investigation_plan.md's "Q3 results" for the real-data cases this
 * was found on, and PHYLIP's protdist.c (makedists()) for the
 * reference approach this follows.
 *
 * Computed via Eigen::Matrix<float,20,20> (Qdecomp.at_eigen_f()/
 * Q_eigen_f()/Q2_eigen_f()) - float, not double: a further ~2x
 * measured on top of the double-Eigen path, error 3-5 orders of
 * magnitude inside CONVERGENCE_TOL (planning/fastprot_ml_speedup_
 * implementation_plan.md). The final reduction below still
 * accumulates in double (cast before summing) - the per-element
 * arithmetic is where float's speed comes from; summing 400 terms is
 * cheap either way and double accumulation costs nothing extra while
 * avoiding adding its own rounding error on top.
 *
 * N is a plain Eigen array here, not a Matrix - likelihood_calc()
 * converts it once, before the Newton loop starts, since N doesn't
 * change across iterations for a given pair and re-reading it via
 * Matrix's bounds-checked accessor on every iteration was pure
 * overhead.
 * @param N The replacement count matrix (N(i,j) = number of observed
 *            replacements from amino acid i to amino acid j),
 *            pre-converted to float Eigen form by likelihood_calc()
 * @param Qdecomp Cached eigendecomposition of Q (see Matrix.hpp)
 * @param t The distance
 * @return The log-likelihood's slope and curvature at t
 */
LikelihoodDerivatives likelihood_slope_curv(const Eigen::Matrix<float, 20, 20> &N, const MatrixExpm &Qdecomp,
                                             double t)
{
    float tf = static_cast<float>(t);
    Eigen::Matrix<float, 20, 20> pt = Qdecomp.at_eigen_f(tf);
    Eigen::Matrix<float, 20, 20> p1t = pt * Qdecomp.Q_eigen_f();
    Eigen::Matrix<float, 20, 20> p2t = pt * Qdecomp.Q2_eigen_f();

    // P(t) is a transition-probability matrix for a rate matrix with
    // strictly positive off-diagonal rates (every model this is used
    // with), so p > 0 for all t > 0 mathematically; floor it
    // defensively against floating-point underflow at extreme t
    // rather than letting a division produce inf/nan. Array
    // expressions (not a hand-written loop with a zero-skip branch,
    // as this used to be) - every N(i,j)==0 cell contributes an exact
    // 0 to the sums regardless (0 * anything), so skipping them saved
    // a division at the cost of a branch; letting Eigen vectorize the
    // whole 20x20 array uniformly is simpler and was measured to not
    // need that branch to be fast.
    Eigen::Array<float, 20, 20> p = pt.array().max(std::numeric_limits<float>::min());
    Eigen::Array<float, 20, 20> ratio = p1t.array() / p;
    Eigen::Array<double, 20, 20> slope_terms = (N.array() * ratio).cast<double>();
    Eigen::Array<double, 20, 20> curv_terms = (N.array() * (p2t.array() / p - ratio.square())).cast<double>();

    return {slope_terms.sum(), curv_terms.sum()};
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
 * @param N The replacement count matrix, already in the fixed 20x20
 *            Eigen form likelihood_slope_curv() needs - callers build
 *            it directly in this form (see ProtDistCalc.cpp's
 *            tally_to_eigen()), so no Matrix round-trip happens
 *            anywhere in the ML path for N.
 * @return The distance with the maximum likelihood
 */
double likelihood_calc(const Eigen::Matrix<float, 20, 20> &N, const MatrixExpm &Qdecomp)
{
    // Eigen's own reductions (unchecked, no bounds-checking overhead)
    // - computed once and reused for both the early-exit check here
    // and kimura_distance()'s starting value, rather than each
    // recomputing them from N separately (a former version did).
    double n_sum = N.sum();
    double n_sum_diag = N.trace();
    if (n_sum - n_sum_diag < DBL_EPSILON)
    {
        return 0;
    }

    double t = kimura_distance(n_sum, n_sum_diag);
    if (t <= 0)
    {
        t = 1;
    }
    double delta = t / 2.0;

    for (int i = 0; i < MAX_ITERATIONS; i++)
    {
        LikelihoodDerivatives d = likelihood_slope_curv(N, Qdecomp, t);

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
