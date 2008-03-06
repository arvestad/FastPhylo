//--------------------------------------------------
//                                        
// File: AML_LeafLifting.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: AML_LeafLifting.hpp,v 1.1 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef AML_LEAFLIFTING_HPP
#define AML_LEAFLIFTING_HPP

#include "SequenceTree.hpp"
#include "dna_pairwise_sequence_likelihood.hpp"

// OPTIMAL LEAF LIFTING
// Returns the likelihood of the optimal leaf lifting.
// The internal nodes of t are filled in with the optimal leaf
// sequences and the edges are given the optimal lengths.
double
computeOptimal_AML_LeafLifting(SequenceTree &t, sequence_model model);

#endif // AML_LEAFLIFTING_HPP
