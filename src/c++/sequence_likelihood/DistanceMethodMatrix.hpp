//--------------------------------------------------
//                                        
// File: DistanceMethodMatrix.hpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: DistanceMethodMatrix.hpp,v 1.4 2006/02/07 14:38:19 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef DISTANCEMETHODMATRIX_HPP
#define DISTANCEMETHODMATRIX_HPP


#include <string>
#include <iostream>
#include "string_compare.hpp"
#include "dna_pairwise_sequence_likelihood.hpp"

//WARNING when we compute the liklihood we don't allow for more than
//prob 1/2 but since these values are going to be used together with
// NJ we need big values for distant permutations
//
// Takes a vector of sequences and computes the pairwise distances.
// The identifiers in the matrix are left untouched.

template <class TreeNode_type>
void
fillDistanceMethodMatrix_P_DISTANCE(std::vector<std::string *> &strvec,
                                    DistanceMatrix< TreeNode_type *, double, Data_init<TreeNode_type*>, Data_printOn<TreeNode_type *>,
                                    Data_init<double>, Data_printOn<double> > &lm){

  lm.resize(strvec.size());
  
  
  for ( size_t i = 0 ; i < strvec.size() ; i++ ){
    int slen =  strvec[i]->size();
    lm.setDistance(i,i, 0);
      
    for ( size_t j = i+1 ; j < strvec.size() ; j++ ){
      TN_string_distance tndist = TN_string_compare(*strvec[i],*strvec[j]);
      simple_string_distance sd = convert_TN_string_distance_to_simple(tndist);
      
      double p_dist = compute_p_distance(slen,sd);
      lm.setDistance(i,j,p_dist);
    }
  }
}


template <class TreeNode_type>
void
fillDistanceMethodMatrix_JC(std::vector<std::string *> &strvec,
                            DistanceMatrix< TreeNode_type *, double, Data_init<TreeNode_type*>, Data_printOn<TreeNode_type *>,
                            Data_init<double>, Data_printOn<double> > &lm){

  lm.resize(strvec.size());
  
  
  for ( size_t i = 0 ; i < strvec.size() ; i++ ){
    int slen =  strvec[i]->size();
    lm.setDistance(i,i, 0);
      
    for ( size_t j = i+1 ; j < strvec.size() ; j++ ){
      TN_string_distance tndist = TN_string_compare(*strvec[i],*strvec[j]);
      simple_string_distance sd = convert_TN_string_distance_to_simple(tndist);
      double p_dist = (compute_JC(slen,sd)).distance;
      lm.setDistance(i,j,p_dist);
    }
  }
}


//
// Takes a vector of sequences and computes the pairwise distances.
//
template <class TreeNode_type>
void
fillDistanceMethodMatrix(std::vector<std::string *> &strvec,
                         DistanceMatrix< TreeNode_type *, double, Data_init<TreeNode_type*>, Data_printOn<TreeNode_type *>,
                         Data_init<double>, Data_printOn<double> > &lm,
                         sequence_model model){

  switch ( model ){
  case P_DISTANCE: fillDistanceMethodMatrix_P_DISTANCE(strvec,lm);
    break;
  case JC: fillDistanceMethodMatrix_JC(strvec,lm);
    break;
  default:
    PROG_ERROR("Unimplemented model " << model);
  }
}



#endif //DISTANCEMETHODMATRIX_HPP
