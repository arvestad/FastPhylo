//--------------------------------------------------
//                                        
// File: AML_SpanningTree.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: AML_SpanningTree.hpp,v 1.1 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef AML_SPANNINGTREE_HPP
#define AML_SPANNINGTREE_HPP

#include <string>
#include <vector>
#include "SequenceTree.hpp"


//Creates a maximum likeliy spanning tree of the input sequences
float
AML_SpanningTree(std::vector<std::string *> &sequences,
                 SequenceTree &resultTree);



//For every internal node a child is added with the same sequence as
//the internal node. This is since in general the MST will be used as
//a start phylogeny.
void
addLeafChildToInternalNodes(SequenceTree &resultTree);


//
// Takes a phylogeny and resolves it into a binary tree
//
void
resolveHighDegreeInternalNodes(SequenceTree &resultTree);

#endif // AML_SPANNINGTREE_HPP
