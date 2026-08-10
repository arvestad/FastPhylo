#pragma once

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
double likelihood_calc(const Matrix &N, const Matrix &Q, const MatrixExpm &Qdecomp);
//! Analytic first and second derivatives of the log-likelihood at t -
//! both computed from the same cached decomposition, so this costs
//! only one extra matrix multiply over evaluating the slope alone.
LikelihoodDerivatives likelihood_slope_curv(const Matrix &N, const Matrix &Q, const MatrixExpm &Qdecomp, double t);
