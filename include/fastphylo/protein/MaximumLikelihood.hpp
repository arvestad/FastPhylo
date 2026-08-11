#pragma once

#include <Eigen/Dense>

// Forward Declarations
class MatrixExpm;

//! First and second derivatives of the log-likelihood with respect to
//! the distance t, both computed analytically from the same cached
//! decomposition - see likelihood_slope_curv().
struct LikelihoodDerivatives
{
    double slope; //! d/dt log L(t)
    double curv;  //! d^2/dt^2 log L(t)
};

//! Calculates the kimura distance for a replacement-count matrix whose
//! total (n_sum) and diagonal (n_sum_diag) sums are n_sum/n_sum_diag -
//! takes the two aggregates rather than the matrix itself since
//! likelihood_calc() already needs both for its own early-exit check
//! and would otherwise compute them a second time here for no reason.
double kimura_distance(double n_sum, double n_sum_diag);
//! Calculates the distance with maximum likelihood, via a safeguarded
//! Newton-Raphson search (see MaximumLikelihood.cpp for the algorithm
//! and why it replaced a finite-difference version that returned
//! unverified, sometimes badly wrong answers on real data - see
//! planning/fastprot_ml_speedup_investigation_plan.md's "Q3 results").
//! N is the replacement-count matrix in the same fixed 20x20 Eigen
//! form likelihood_slope_curv() needs - callers build it directly in
//! this form (see ProtDistCalc.cpp's tally_to_eigen()) rather than
//! via the general Matrix class, which N never otherwise needs.
//! Qdecomp must be MatrixExpm(Q) - passed in rather than recomputed
//! here because Q is the same rate matrix for every pair in a
//! calculate_ml_dists() run, and decomposing it is the expensive part
//! (see Matrix.hpp).
double likelihood_calc(const Eigen::Matrix<float, 20, 20> &N, const MatrixExpm &Qdecomp);
//! Analytic first and second derivatives of the log-likelihood at t,
//! computed from the cached decomposition's float-precision copy (see
//! Matrix.hpp's MatrixExpm - planning/fastprot_ml_speedup_
//! implementation_plan.md).
LikelihoodDerivatives likelihood_slope_curv(const Eigen::Matrix<float, 20, 20> &N, const MatrixExpm &Qdecomp,
                                             double t);
