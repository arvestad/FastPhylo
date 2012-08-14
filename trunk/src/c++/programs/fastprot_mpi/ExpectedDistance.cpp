#include "ExpectedDistance.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <cstdio>
#include "Matrix.hpp"

  double op_pow2(double i){return std::pow(i, 2);}
  double op_pow3(double i){return std::pow(i, 3);}
  double op_pow4(double i){return std::pow(i, 4);}
  
  /*!
   * If needed, precomputes a vector with distance samples,
   * vectors with the distance samples ², ³ and ⁴. Calculates the prior 
   * probability and prior distribution.
   * @param tp The type of prior probability
   * @param step_size The step size that should be used between distance
   *         samples
   * @param Q Rate matrix for replacements from one amino acid to another
   * @param eq Equilibrium frequencies for amino acids
   */
  void initialize_ed(const type_prior &tp, int step_size, const Matrix &Q, const DblVec &eq){
    static bool done = false; // We don't want to precompute values several times

    if (!done) {
      // Precomputation of various values

      speed = step_vec[step_size-1];
      DSamples.reserve(static_cast<int>(ceil((max_distance + 1)/speed)));
      // Creation of the vector with distance samples
      for (int i=1; i<=max_distance; i+=speed){
        DSamples.push_back(i);
      }

      // Amount of distance samples, size of a lot of vectors
      nr_distances = DSamples.size(); 
      
      DSamples_prev.resize(nr_distances); 
      std::rotate_copy(DSamples.begin(), DSamples.end()-1, DSamples.end(), DSamples_prev.begin());
      DSamples_prev[0] = 0;
 
      // Calculate SampleGap 
      sample_g1 = log(DSamples[0]/2.0);
      sample_g_other = log(speed/2.0);
      sample_g_last = log( (DSamples.back()-DSamples_prev.back()) / 2.0);

      // Calculate the prior probability
      prior_prob = prior_probability(Q, eq);


      // Calculate normal or flat prior distribution
      switch(tp) {
        case norm: 
          log_prior_dist = log_norm_prior();
          break;
        case flat:
          log_prior_dist = DblVec(nr_distances, 0); //log(1) = 0
          break;
      }
      done = true;
    }
  }

  /*!
   * Function for calculating the expected distance between
   * two protein sequences.
   * @param N The replacement count matrix, N(i,j) contains the number of actual
   *            replacements from amino acid i to amino acid j
   * @return The expected distance between two protein sequences
   */
double calculate_ed(const Matrix &N){
  return expected_distance(N);
}
  
  /*!
   * Calculates the logarithm of the prior probability matrix of a set of
   * replacements at distance d = DSamples[i].
   * The probability is /f$ P(d) = ln( diag(eq) * e^{Q*d} ) /f$ 
   * @param Q Rate matrix for replacements from one amino acid to another
   * @param eq Equilibrium frequencies for amino acids
   * @return A vector with matrices P(d)
   */
MatVec prior_probability(const Matrix &Q, const DblVec &eq){
  // Return vector with matrices
  // Each matrix m = log(diag(eq)*expm(Q*d))
  MatVec pVec;
  pVec.reserve(nr_distances);

  // Matrix with the exponentials of the products Q*DSamples[i]
  MatVec exp_vec = Q.expm(DSamples);

  // Creating a diagonal matrix of eq
  Matrix d_eq(eq);

  // Calculating the prior probability P(d) = ln(d_eq*e^(Q*d))
  for(MatVec::const_iterator it = exp_vec.begin(); it != exp_vec.end(); it++){
    pVec.push_back(Matrix::mult(d_eq, *it));
  }

  for (MatVec::iterator it = pVec.begin(); it != pVec.end(); it++)
    it->mlog();

  return pVec;   

}

  /*! 
   * Calculates the distribution function for the normal distribution
   * with mean value = norm_mean and variance = norm_var
   * @return A vector with values of the distribution function for the normal 
   *          distribution at distance DSamples[i]
   */
DblVec log_norm_prior(){
  const double pi = 3.141592653589793238462643383;
  double x = 1.0/(norm_var*sqrt(2.0*pi));
  double divisor = 2.0*pow(norm_var, 2);

  DblVec result;
  result.reserve(nr_distances);
  for (DblVec::const_iterator it = DSamples.begin(); it != DSamples.end(); it++){
    double dividend = pow((*it - norm_mean), 2);
    result.push_back(log(x*exp(-(dividend/divisor))));
  }
  return result; 
}
  
  /*!
   * Function calculating the expected distance 
   * @param N The replacement count matrix, N(i,j) contains the number of actual
   *            replacements from amino acid i to amino acid j
   * @param sd True if the standard deviation should be calculated, false otherwise
   * @return The expected distance
   */
double expected_distance(const Matrix &N){
  // Check if there are data available
  if (N.sum() == 0){
    return -1;
    std::cout << "No data!" << std::endl;
  }

  // Are the sequences identical? If so, there should only be values on the 
  // main diagonal
  if (N.sum() == N.sum_diag())
    identical = 1;
  else 
    identical = 0;

  // Calculate the posterior probability
  DblVec post_prob = posterior_probability(N);

  return integrate(post_prob);
  
  /*DblVec fnk(nr_distances);

  // fnk[i] = loglikelihood(N*P(i)) + log_prior_distribution(i)
  // PROF: nr_distances = DSample.size() ~ 400/size, where size: 1...10
  
  double term = 0.0;
 #pragma omp parallel default(none) shared(fnk, term, prior_prob, log_prior_dist, speed, nr_distances, identical)
  {
    //printf("Number of threads: %d\n", omp_get_thread_num());
 #pragma omp for
    for (int i=0; i<nr_distances; i++){
      fnk[i] = elem_mult(N, prior_prob[i]).sum() + log_prior_dist[i];
      //    fnk.push_back( (elem_mult(N, prior_prob[i])).sum() + log_prior_dist[i]);
    }

    // Intergration  with trapezoid rule for total posterior probability
    double ptot;
    if (identical)                            // Special treatment for the first value
      ptot = log(1+exp(fnk.front())) - log(2); 
    else
      ptot = fnk.front() - log(2);
    // For the rest of the values
    //#pragma omp for reduction(+: ptot)    
    // XXX: How to take care of log_add addition?
    for (int i=1; i<nr_distances; i++){
      double term = log_add(fnk[i], fnk[i-1]) + sample_gap(i);
      ptot = log_add(ptot, term);
    }

    DblVec res(nr_distances, 0);
#pragma omp for
    for (int i=0; i<nr_distances; i++){
      res[i] = exp( fnk[i] - ptot);
    }

    // Integrate 
    if (identical)                  // Special treatment for the first value
      term = (1.0+res.front()); 
    else
      term = res.front();
#pragma omp for reduction(+:term)
    for (int i=2; i<nr_distances; i++)
      term += (res[i]*i + res[i-1]*(i-1)) * speed*speed;
  } // End of parallel region
  term = term/2.0;
  return term;*/
}

  /*!
   * Compute the posterior probability of observing a set of replacements 
   * The integration code demands that the samples are uniformly distributed.
   * Numerical integration using simple linear interpolation. 
   * @param N The replacement count matrix, N(i,j) contains the number of actual
   *            replacements from amino acid i to amino acid j
   * @return A vector with the values of the posterior probability of observing
   *          a set of replacements
   */
DblVec posterior_probability(const Matrix &N){

  DblVec fnk(nr_distances);

  // fnk[i] = loglikelihood(N*P(i)) + log_prior_distribution(i)
  // PROF: nr_distances = DSample.size() ~ 400/size, where size: 1...10

//#pragma omp parallel for default(none) private(i,nr_distances) shared(fnk, N, prior_prob, log_prior_dist)
  for (int i=0; i<nr_distances; i++){
    fnk[i] = elem_mult(N, prior_prob[i]).sum() + log_prior_dist[i];
    //    fnk.push_back( (elem_mult(N, prior_prob[i])).sum() + log_prior_dist[i]);
  }

  // Intergration  with trapezoid rule for total posterior probability
  double ptot;
  if (identical)                            // Special treatment for the first value
    ptot = log(1+exp(fnk.front())) - log(2); 
  else
    ptot = fnk.front() - log(2);
  // For the rest of the values
  for (int i=1; i<nr_distances; i++){
    double term = log_add(fnk[i], fnk[i-1]) + sample_gap(i);
    ptot = log_add(ptot, term);
  }

  DblVec res(nr_distances);
  for (int i=0; i<nr_distances; i++) {
      res[i] = exp(fnk[i] - ptot);
  }

  return res;

}

/*
 * Integration with the trapezoid rule 
 * @param fnk Vector containing the function to be integrated
 * @return The value of the integral
 */
double integrate(const DblVec fnk){
  double term = 0.0;
  if (identical)                  // Special treatment for the first value
    term = (1.0+fnk.front()); 
  else
    term = fnk.front();
  for (int i=2; i<nr_distances; i++)
    term += (fnk[i]*i + fnk[i-1]*(i-1)) * speed*speed;
  return term/2.0;
}
