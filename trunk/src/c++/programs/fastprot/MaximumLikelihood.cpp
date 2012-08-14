#include "MaximumLikelihood.hpp"
#include <cmath>
#include <cfloat>
#include <iostream>
#include "Matrix.hpp"

  /*! 
   * Kimura distance, not corrected, used as a starting value for Newton's
   * method in likelihood_calc()
   * @param N The replacement count matrix, N(i,j) contains the number of actual
   *            replacements from amino acid i to amino acid j
   * @return The kimura distance in PAMs for the matrix N
   */
  double kimura_distance(const Matrix &N){
    double n_sum = N.sum();
    double d = (n_sum-N.sum_diag())/n_sum;

    double adjusted = d + 0.2*pow(d, 2);
    if (adjusted > 0.854) // Infinite distance
      adjusted = 0.854;

    adjusted = -100 * log(1-adjusted);
    return adjusted;
  }

  /*! 
   * Computes the distance with the maximum likelihood using 
   * Newton-Rhapson
   * @param N The replacement count matrix, N(i,j) contains the number of actual
   *            replacements from amino acid i to amino acid j
   * @param Q Rate matrix for replacements from one amino acid to another
   * @return The distance with the maximum likelihood
   */
  double likelihood_calc(const Matrix &N, const Matrix &Q){
    if (N.sum() - N.sum_diag() < DBL_EPSILON)
      return 0;
    
    // Starting value for t
    double t = kimura_distance(N);
    
    if (t == 0)
      t = 1;

    // Newton-Rhapson
    double delta, tol;
    delta = tol = 0.001;
    int maxit = 50; // max iterations

    double l_d = likelihood_deriv(N, Q, t);
    for (int i=0; i<maxit; i++){
      if (fabs(l_d) < tol){ // If the derivative is small enough
        return t;
      }
      double l_new = likelihood_deriv(N, Q, t+delta);
      double deriv = (l_new - l_d) / delta;
      t = t - l_d / deriv;
      if (t < 1) // 1 is the smallest possible distance
        return 1;
      if (t > 500) // 500 is infinity
        return 500;
      if (fabs(l_d) < fabs(l_new)) // The derivative is getting larger..
        return t;
      l_d = l_new;
    }
    return t;

  }
  /*!
   * Derivative of the loglikelihood
   * log L(t) = log \prod_{i,j} p_{i,j}^{N_{i,j}}(t) = \sum N_{i,j} log(p_{i,j}(t))
   * (log L(t))' = \sum_{i,j} N_{i,j} \cdot \frac{1}{p_{i,j}(t)} \cdot p_{i,j}'t()
   * Which means that it becomes
   * log L(t) = \sum N .* P'(t) ./ P(t) 
   * where
   * P(t) = e^{Qt} 
   * P'(t) = Q e^{Qt} 
   * @param N The replacement count matrix, N(i,j) contains the number of actual
   *            replacements from amino acid i to amino acid j
   * @param Q Rate matrix for replacements from one amino acid to another
   * @param t The distance
   * @return The derivative
   */
  double likelihood_deriv(const Matrix &N, const Matrix &Q, double t){

    // P(t) = e^Qt
    Matrix pt = Q.expm(DblVec(1, t))[0];
    
    // P'(t) = e^Qt * (Qt)' 
    Matrix direction = Matrix::mult(pt, Q);
    
    Matrix temp = elem_mult(N, direction);
    Matrix temp2 = elem_div(temp, pt);

    return temp2.sum();
  }
