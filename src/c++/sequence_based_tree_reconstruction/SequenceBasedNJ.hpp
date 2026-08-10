//--------------------------------------------------
//
// File: SequenceBasedNJ.hpp
//
// Author: Isaac Elias
// e-mail: isaac@nada.kth.se
//
//--------------------------------------------------
#ifndef SEQUENCEBASEDNJ_HPP
#define SEQUENCEBASEDNJ_HPP

#include <vector>
#include "fastphylo/core/SequenceTree.hpp"

void computeSequenceBasedNJ(std::vector<Sequence> &seqs, SequenceTree &resultTree);

#endif // SEQUENCEBASEDNJ_HPP
