#ifndef _MAXLIKE_HPP_
#define _MAXLIKE_HPP_

  // Forward Declarations
  class Matrix;
  class MatrixExpm;

  //! Calculates the kimura distance of matrix N
  double kimura_distance(const Matrix &N);
  //! Calculates the distance with maximum likelihood. Qdecomp must be
  //! MatrixExpm(Q) - passed in rather than recomputed here because Q is
  //! the same rate matrix for every pair in a calculate_ml_dists() run,
  //! and decomposing it is the expensive part (see Matrix.hpp).
  double likelihood_calc(const Matrix &N, const Matrix &Q, const MatrixExpm &Qdecomp);
  //! Calculates the derivative of likelihood
  double likelihood_deriv(const Matrix &N, const Matrix &Q, const MatrixExpm &Qdecomp, double t);

#endif
