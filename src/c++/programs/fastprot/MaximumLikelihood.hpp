#ifndef _MAXLIKE_HPP_
#define _MAXLIKE_HPP_

  // Forward Declarations
  class Matrix;
  
  //! Calculates the kimura distance of matrix N
  double kimura_distance(const Matrix &N);
  //! Calculates the distance with maximum likelihood
  double likelihood_calc(const Matrix &N, const Matrix &Q);
  //! Calculates the derivative of likelihood
  double likelihood_deriv(const Matrix &N, const Matrix &Q, double t);

#endif
