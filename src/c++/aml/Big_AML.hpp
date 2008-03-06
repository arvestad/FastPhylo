//--------------------------------------------------
//                                        
// File: Big_AML.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: Big_AML.hpp,v 1.1 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef BIG_AML_HPP
#define BIG_AML_HPP

#include "dna_pairwise_sequence_likelihood.hpp"
#include "SequenceTree.hpp"
#include <vector>
#include "Sequence.hpp"
#include <string>
#include "InitAndPrintOn_utils.hpp"
#include "DistanceMatrix.hpp"


//
// Takes as input a vector of sequences and computes an AML tree using
// the spanning tree and NJ algorithm. This tree has a guaranteed
// 2-approx ratio and is also guaranteed to return the correct tree
// for nearly additive data.  The ancestoral sequences are first given
// by the MST and then improved through local search.

void
Big_AML(std::vector<Sequence> &seqs, SequenceTree &resultTree, sequence_model m);


//
// Creates a maximum likeliy spanning tree of the sequences in the
// identifiers.  To compute a regular MST tree the each input tree
// node should belong to its own unique tree, these trees will then be
// joined into one tree single tree.  If the nodes don't belong to
// unique trees a spanning tree will be created between the nodes as
// long as the nodes are in seperate trees.
//
// The sum of loglikelihoods of the inferred edges is returned.
//
double
AML_SpanningTree(SequenceTree::NodeMatrix &lm);

//
// For each node in the vector with degree > 1 a leaf child is added with 
// the same sequence data as the internal node.
//
void
addLeafChildToInternalNodesInVector(std::vector<SequenceTree::Node *> &nodes);


#endif // BIG_AML_HPP
