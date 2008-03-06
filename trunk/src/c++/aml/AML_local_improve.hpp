//--------------------------------------------------
//                                        
// File: AML_local_improve.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: AML_local_improve.hpp,v 1.1 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef AML_LOCAL_IMPROVE_HPP
#define AML_LOCAL_IMPROVE_HPP

#include "SequenceTree.hpp"
#include "dna_pairwise_sequence_likelihood.hpp"

void
AML_local_improve(SequenceTree &t, sequence_model m);

#endif // AML_LOCAL_IMPROVE_HPP
