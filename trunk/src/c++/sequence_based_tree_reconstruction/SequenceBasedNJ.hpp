//--------------------------------------------------
//                                        
// File: SequenceBasedNJ.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: SequenceBasedNJ.hpp,v 1.1 2006/03/31 17:11:33 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef SEQUENCEBASEDNJ_HPP
#define SEQUENCEBASEDNJ_HPP

#include <vector>
#include "SequenceTree.hpp"

void
computeSequenceBasedNJ(std::vector<Sequence> &seqs, SequenceTree &resultTree);

#endif // SEQUENCEBASEDNJ_HPP
