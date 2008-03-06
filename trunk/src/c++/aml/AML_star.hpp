//--------------------------------------------------
//                                        
// File: AML_star.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: AML_star.hpp,v 1.1 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef AML_STAR_HPP
#define AML_STAR_HPP

#include <string>
#include "dna_pairwise_sequence_likelihood.hpp"

//Fills in center with the most likely center string
//of a b and c. Returns the likelihood of the star
double find_AML_star(std::string &a, 
	       std::string &b, 
	       std::string &c,
	       std::string &center,
               sequence_model model);

//
//computes the likelihood of the given strings
// 
bool
improve_AML_star(std::string &a, 
                 std::string &b, 
                 std::string &c,
                 std::string &center,
                 sequence_model model);

#endif // AML_STAR_HPP
