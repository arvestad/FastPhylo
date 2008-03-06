//--------------------------------------------------
//                                        
// File: string_compare.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: string_compare.hpp,v 1.1 2006/01/25 08:03:37 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef STRING_COMPARE_HPP
#define STRING_COMPARE_HPP


#include <string>
#include "dna_pairwise_sequence_likelihood.hpp"
#include <iostream>


//Straigt forward hamming distance. Doesn't check for ambiguities or
//gaps.
int hamming_distance(const std::string &s1,
		     const std::string &s2);

//Returns sufficient statistica for the Tamura Nei model
//i.e. the number of purine transitions, pyrimidine transitions, transversions, and the number of deletions
TN_string_distance
TN_string_compare(const std::string &s1,
                  const std::string &s2);

//returns the number of deleted characters.
int 
complete_dna_string_compare(float divergence_matrix[4][4],//a 4x4 matrix will be filled in with the frequences
                            const std::string &s1,
                            const std::string &s2);


TN_string_distance
divergence_matrix_2_TN_distance(float divergence_matrix[4][4], int deleted);


void
print_divergence_matrix(std::ostream &out, float div_matrix[4][4]);


#endif // STRING_COMPARE_HPP











