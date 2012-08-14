#ifndef _EXPECTMAX_HPP_
#define _EXPECTMAX_HPP_

#include <vector>
#include <cmath>
#include "Matrix.hpp"

  // Forward Declarations
  class Matrix;

  typedef std::vector<double> DblVec;
  typedef std::vector<Matrix> MatVec;
  
  //! Enum for the type of prior probability distribution
  //! Norm - normal distribution, flat - flat distribution
  enum type_prior{norm, flat};

      //! Precomputes a number of values if needed.
      void initialize_ed(const type_prior &tp, int step_size, 
          const Matrix &Q, const DblVec &eq);
      //! Calculates the ED with standard deviation
      double calculate_ed_with_sd(const Matrix &N);
      //! Calculates the ED without standard deviation
      double calculate_ed(const Matrix &N);
      //! Access method for the standard deviation
      double get_standard_deviation();

  
      // Computational constants

      //! Maximum distance
      static const int max_distance = 400;
      //! Minimum distance
      static const int min_distance = 1;
      // This describes the prior distribution for, for example, 
      // pairwise distances in Pfam
      //! Mean value for the normal distribution
      static const int norm_mean = 115;
      //! Variance for the normal distribution
      static const int norm_var = 50;
      //! Array for storing different step sizes
      static const int step_vec[] = {1,2,4,5,8,10,20,25};
      //static const int step_vec[] = {4,10,25,40,50,80,100,200};
      //static const int step_vec[] = {5,10,20,25,40,80,100,200};
      
      // Pre-calculated values

      //! Vector with calculated values for the prior distribution
      static DblVec log_prior_dist;
      //! Vector with matrices containing the prior probability (of observing
      //! a set of replacements) 
      static MatVec prior_prob; 
      //! Vectors with distance samples
      static DblVec DSamples, DSamples2, DSamples3, DSamples4;
      //! Vectors with shifted values of distance samples
      static DblVec DSamples_prev, DSamples2_prev, DSamples3_prev, DSamples4_prev;
      //! The number of distance samples (DSamples.size())
      static int nr_distances;
      //! Values for sample gaps 
      static double sample_g1;
      static double sample_g_other;
      static double sample_g_last;
      //! The step size between distance samples
      static double speed;


      //! Is the two sequences identical?
      static int identical;
      //! Standard deviation for the ED-value
      static double std_deviation;
      

      // Help functions
      //! Calculates the prior probability matrix
      MatVec prior_probability(const Matrix &Q, const DblVec &eq);
      //! Calculates the distribution function for the normal distribution
      DblVec log_norm_prior();
      //! Calculates the posterior probability
      DblVec posterior_probability(const Matrix &N);
      
      //! Calculates the expected distance 
      double expected_distance(const Matrix &N, bool sd);
      
      //! Returns the logarithm of the sum e^lhv and e^rhv
      double log_add(const double lhv, const double rhv);
      //! Returns the sample gap between DSample[i] and DSample_prev[i]
      double sample_gap(int i);
      //! Integrates fnk over values in DSamples
      double integrate(const DblVec fnk);


  /*!
   * Returns the logarithm of the gap between two distance samples.
   * Makes use of the fact that since a constant step size is used
   * the sample gap between most of the distances is the same.
   * As the code is written right now sample_gap(0) shouldn't be used,
   * therefore it is placed last.
   * This function can be replaced by a vector instead, which contains
   * all the sample gaps, but it would contain a lot of duplicate values
   * (all values except for the last and the first actually).
   * @param i The index
   * @return The logarithm of the half sample gap between DSample[i] and
   *            DSample_prev[i]
   */
  inline double sample_gap(int i){
    if (i == nr_distances-1)
      return sample_g_last;
    else if (i != 0)
      return sample_g_other;
    else
      return sample_g1;
  }

  /*! 
   * Adds two values. 
   * @param lhv Logarithm of value x
   * @param rhv Logarithm of value y
   * @return The logarithm of the sum x*y
   */
  inline double log_add(const double lhv, const double rhv){
    double result;
    if (lhv < rhv) 
      result = rhv + log(1 + exp(lhv-rhv));
    else
      result = lhv + log(1 + exp(rhv-lhv));

    return result;
  }

#endif
