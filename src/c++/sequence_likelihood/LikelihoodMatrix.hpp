//--------------------------------------------------
//                                        
// File: LikelihoodMatrix.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: LikelihoodMatrix.hpp,v 1.9 2006/12/08 11:09:14 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef LIKELIHOODMATRIX_HPP
#define LIKELIHOODMATRIX_HPP

#include "DistanceMatrix.hpp"
#include <string>
#include <vector>
#include "InitAndPrintOn_utils.hpp"
#include "dna_pairwise_sequence_likelihood.hpp"
#include <iostream>
#include "string_compare.hpp"

//
// Takes a vector of sequences and computes the pairwise liklihoods.
// The identifiers in the matrix are left untouched.
// The likelihood is the same under Jukes Cantor
template <class TreeNode_type>
void
fillLikelihoodMatrix_P_DISTANCE(std::vector<std::string *> &strvec,
                                DistanceMatrix< TreeNode_type *, double, Data_init<TreeNode_type*>, Data_printOn<TreeNode_type *>,
                                Data_init<double>, Data_printOn<double> > &lm){

  lm.resize(strvec.size());
  
  double hamdist;
  
  
  for ( size_t i = 0 ; i < strvec.size() ; i++ ){
    double slen = 1.0* strvec[i]->size();
    lm.setDistance(i,i, 0);
      
    for ( size_t j = i+1 ; j < strvec.size() ; j++ ){
      hamdist = hamming_distance(*(strvec[i]),*(strvec[j]));
      double optprob = hamdist/slen;
      optprob = ( 0.499 < optprob ? 0.499 : optprob );//we don't allow for more than 0.5 probability of change
      if ( hamdist != 0 )
        lm.setDistance(i,j, hamdist*log(optprob/3.0) + (slen - hamdist)*log(1 - optprob));
      else
        lm.setDistance(i,j,0.0);

    }
  }
}


//
// Takes a vector of sequences and computes the pairwise liklihoods.
//
template <class TreeNode_type>
void
fillLikelihoodMatrix(std::vector<std::string *> &strvec,
                     DistanceMatrix< TreeNode_type *, double, Data_init<TreeNode_type*>, Data_printOn<TreeNode_type *>,
                     Data_init<double>, Data_printOn<double> > &lm,
                     sequence_model model){

  switch ( model ){
    //likelihood for p-distance and JC is the same
  case P_DISTANCE: case JC:fillLikelihoodMatrix_P_DISTANCE(strvec,lm);
    break;
  default:
    PROG_ERROR("Unimplemented model " << model);
  }
}

//
// FIX FACTOR CORRECTION
//
// All NaNs are set to minValInMatrix * factor
// Returns true if some value was corrected
template <class T>
bool
fixFactorCorrection_ToMIN(DistanceMatrix< T, double, Data_init<T>, Data_printOn<T>,
                        Data_init<double>, Data_printOn<double> > &lm,
                        double fixFactor = 2.0){
  
  double minVal = lm.getDistace(0,0);
  bool shouldCorrect = false;
  
  for ( size_t i = 0 ; i < lm.getSize() ; i++ ){
      
    for ( size_t j = i+1 ; j < lm.getSize() ; j++ ){
      double dist = lm.getDistance(i,j);
      if ( isnan(dist) )
        shouldCorrect = true;
      else{
        if ( dist < minVal )
          minVal = dist;
      }
    }
  }

  if ( ! shouldCorrect )
    return false;

  minVal = minVal * fixFactor;
  for ( size_t i = 0 ; i < lm.getSize() ; i++ ){
      
    for ( size_t j = i+1 ; j < lm.getSize() ; j++ ){
      double dist = lm.getDistance(i,j);
      if ( isnan(dist) ){
        PRINT_V(minVal);
        lm.setDistance(i,j,minVal);
      }
    }
  }

  return true;
}

template <class T>
bool
fixFactorCorrection_ToMAX(DistanceMatrix< T, double, Data_init<T>, Data_printOn<T>,
                        Data_init<double>, Data_printOn<double> > &lm,
                        double fixFactor = 2.0){
  
  double maxVal = lm.getDistance(0,0);
  bool shouldCorrect = false;
  
  for ( size_t i = 0 ; i < lm.getSize() ; i++ ){
      
    for ( size_t j = i+1 ; j < lm.getSize() ; j++ ){
      double dist = lm.getDistance(i,j);
      if ( isnan(dist) )
        shouldCorrect = true;
      else{
        if ( dist > maxVal )
          maxVal = dist;
      }
    }
  }

  if ( ! shouldCorrect )
    return false;

  maxVal = maxVal * fixFactor;
  for ( size_t i = 0 ; i < lm.getSize() ; i++ ){
      
    for ( size_t j = i+1 ; j < lm.getSize() ; j++ ){
      double dist = lm.getDistance(i,j);
      if ( isnan(dist) ){
        lm.setDistance(i,j,maxVal);
      }
    }
  }

  return true;
}

#endif // LIKELIHOODMATRIX_HPP




















