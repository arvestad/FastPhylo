//--------------------------------------------------
//                                        
// File: AML_given_edge_probabilities.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: AML_given_edge_probabilities.hpp,v 1.1 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef AML_GIVEN_EDGE_PROBABILITIES_HPP
#define AML_GIVEN_EDGE_PROBABILITIES_HPP

#include "SequenceTree.hpp"
#include "dna_pairwise_sequence_likelihood.hpp"


double
computeAML_given_edge_probabilities(SequenceTree &t);


void
convert_edge_lengths_to_probabilities(SequenceTree &t, sequence_model model);


#endif // AML_GIVEN_EDGES_HPP
