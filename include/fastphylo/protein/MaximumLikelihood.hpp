#pragma once

#include <Eigen/Dense>

// Forward Declarations
class Matrix;
class MatrixExpm;

//! First and second derivatives of the log-likelihood with respect to
//! the distance t, both computed analytically from the same cached
//! decomposition - see likelihood_slope_curv().
struct LikelihoodDerivatives
{
    double slope; //! d/dt log L(t)
    double curv;  //! d^2/dt^2 log L(t)
};

//! Calculates the kimura distance of matrix N
double kimura_distance(const Matrix &N);
//! Calculates the distance with maximum likelihood, via a safeguarded
//! Newton-Raphson search (see MaximumLikelihood.cpp for the algorithm
//! and why it replaced a finite-difference version that returned
//! unverified, sometimes badly wrong answers on real data - see
//! planning/fastprot_ml_speedup_investigation_plan.md's "Q3 results").
//! Qdecomp must be MatrixExpm(Q) - passed in rather than recomputed
//! here because Q is the same rate matrix for every pair in a
//! calculate_ml_dists() run, and decomposing it is the expensive part
//! (see Matrix.hpp).
double likelihood_calc(const Matrix &N, const MatrixExpm &Qdecomp);
//! Analytic first and second derivatives of the log-likelihood at t,
//! computed from the cached decomposition's float-precision copy (see
//! Matrix.hpp's MatrixExpm - planning/fastprot_ml_speedup_
//! implementation_plan.md) - N must already be converted to the same
//! fixed 20x20 form (likelihood_calc() does this once per pair, not
//! once per call - N doesn't change across Newton iterations).
LikelihoodDerivatives likelihood_slope_curv(const Eigen::Matrix<float, 20, 20> &N, const MatrixExpm &Qdecomp,
                                             double t);
