//--------------------------------------------------
//                                        
// File: SequenceTree_MostParsimonious.cpp                              
//         
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: SequenceTree_MostParsimonious.cpp,v 1.2 2006/12/19 09:24:36 isaac Exp $                                 
//
// IMPLEMENTATION of SequenceTree::computeMostParsimoniousSequences()
//
//--------------------------------------------------

#include "fastphylo/core/SequenceTree.hpp"
#include <string>
#include "fastphylo/core/nucleotide.hpp"
#include <vector>
#include "fastphylo/core/log_utils.hpp"
#include "fastphylo/core/std_c_utils.h"

using namespace std;


//The parsimony information is saved in a vector of 4 integers. Each position
//has one parsimony inforamtion
using p_info = vector<size_t>;

using p_vector = vector<p_info>;


std::istream &
operator>>(std::istream &in,p_info p){
  in >>p[DNA_A_]; in >> p[DNA_C_];in >> p[DNA_G_];in >> p[DNA_T_];
  return in;
}

std::ostream&
operator<<(std::ostream & os, p_info  p){
  os << p[DNA_A_]; os << p[DNA_C_]; os << p[DNA_G_]; os << p[DNA_T_];
  return os;
}

using ParsimonyTree = Tree<p_vector>;

static void
backtrack(SequenceTree::Node *sn, ParsimonyTree::Node *pn);

namespace {
// Fitch-parsimony bottom-up DP core: for node's data at position j,
// computes the cost of assigning each of the 4 symbols to this
// position, given its children's already-computed costs (min over
// each child of that child's cost for a symbol, plus 1 if it differs
// from the symbol being priced here). The triple-nested loop
// (symbol x child x child-symbol) below was the main contributor to
// computeMostParsimoniousSequences()'s cognitive complexity (45).
void computeParsimonyAtPosition(ParsimonyTree::Node *node, size_t j) {
  constexpr size_t numSymbols = 4;
  p_info &p = node->data[j];
  for ( size_t sym = 0 ; sym < numSymbols ; sym++ ){
    p[sym] = 0;
    ParsimonyTree::Node *child = node->getRightMostChild();
    for ( ; child != nullptr ; child = child->getLeftSibling() ){
      p_info &child_p = (child->data)[j];
      size_t score = child_p[0] + (sym != 0 ? 1 : 0);
      for ( size_t symC = 1 ; symC < numSymbols ; symC++ ){
        size_t symCp = child_p[symC] + (sym != symC ? 1 : 0);
        score = MIN(score, symCp);
      }
      //add the best cost the the parsimony at the parent
      p[sym] += score;
    }//end children loop
  }//end symbol loop
}
} // namespace

size_t
SequenceTree::computeMostParsimoniousSequences(){

  SequenceTree::Node *oldRoot = getRoot();
  reRootAtHighDegreeNode();

  p_vector defaultVec;
  ParsimonyTree pTree(*this,defaultVec);
  
  
  ParsimonyTree::NodeVector nodes;
  nodes.reserve(pTree.getNumNodes());

  //init the leafs
  p_info defaultVal(4);
  defaultVal[0]=10;defaultVal[1]=10;defaultVal[2]=10;defaultVal[3]=10;
  SequenceTree::NodeVector origLeafs;
  origLeafs.reserve(getNumLeafs());
  addLeafs(origLeafs);
  pTree.addLeafs(nodes);
  string &seq = oldRoot->data.s.seq;
  for ( size_t i = 0 ;i < origLeafs.size() ; i++ ){
    seq = origLeafs[i]->data.s.seq;
    p_vector &pos = nodes[i]->data;
   
    pos.resize(seq.length(),defaultVal);
   
    for ( size_t j = 0 ; j < seq.length() ; j++ ){
      nucleotide n = char2nucleotide(seq[j]);
      pos[j][n] = 0;
    }
  }

  //constants
  const size_t seqlen = seq.length();
  const size_t numNodes = getNumNodes();


  //BOTTOM UP
  nodes.clear();
  pTree.addNodesInPostfixOrder(nodes);

  for ( size_t i = 0 ; i < numNodes ; i++ ){
    if ( nodes[i]->isLeaf() ) { continue;
}
    p_vector &pos = nodes[i]->data;
    pos.resize(seq.length(),defaultVal);
    //for each position in the sequence, compute the parsimony score
    //for each possible symbol by taking the minimum sum over the children
    for( size_t j = 0 ; j < seqlen ; j++ ){
      computeParsimonyAtPosition(nodes[i], j);
    }//end sequence position
  }//end node loop


  //BACK TRACK
  //first create the best sequence in the root and then backtrack down.
  SequenceTree::Node *sroot = getRoot();
  ParsimonyTree::Node *proot = pTree.getRoot();
  p_vector &pos = proot->data;
  seq = sroot->data.s.seq;
  seq.resize(seqlen,'x');
  size_t totalParsimony = 0;
  for ( size_t i = 0 ; i < seqlen ; i++ ){
    p_info &p = pos[i];
    nucleotide min_n = DNA_A_;
    size_t score = p[DNA_A_];
    if ( p[DNA_C_] < score ){
      min_n = DNA_C_;score = p[DNA_C_];
    }
    if ( p[DNA_G_] < score ){
      min_n = DNA_G_;score = p[DNA_G_];
    }
    if ( p[DNA_T_] < score ){
      min_n = DNA_T_;score = p[DNA_T_];
    }
    seq[i] = nucleotide2char(min_n);
    totalParsimony+= score;
  }
  SequenceTree::Node *child = sroot->getRightMostChild();
  ParsimonyTree::Node *childP = proot->getRightMostChild();
  for ( ; child != nullptr ; child = child->getLeftSibling(), childP = childP->getLeftSibling() ) {
    backtrack(child,childP);
}
  
  //reroot the tree to its old position
  reRootAt(oldRoot);

  return totalParsimony;
}



static void
backtrack(SequenceTree::Node *sn, ParsimonyTree::Node *pn){

  if ( sn->isLeaf() ) {
    return;
}

  
  p_vector &pos = pn->data;
  const size_t seqlen = pos.size();
  string &seq = sn->data.s.seq;
  seq.resize(seqlen,'x');
  string &parentseq = sn->getParent()->data.s.seq;

  for ( size_t i = 0 ; i < seqlen ; i++ ){
    p_info &p = pos[i];
    nucleotide parent_n = char2nucleotide(parentseq[i]);
    nucleotide min_n = parent_n;
    size_t score = p[min_n];
    if ( p[DNA_A_]+1 < score ){
      min_n = DNA_A_;
      score = p[DNA_A_]+1;
    }
    if ( p[DNA_C_]+1 < score ){
      min_n = DNA_C_;
      score = p[DNA_C_]+1;
    }
    if ( p[DNA_G_]+1 < score ){
      min_n = DNA_G_;
      score = p[DNA_G_]+1;
    }
    if ( p[DNA_T_]+1 < score ){
      min_n = DNA_T_;
    }
    seq[i] = nucleotide2char(min_n);
  }

  

  //continue backtracking
  SequenceTree::Node *child = sn->getRightMostChild();
  ParsimonyTree::Node *childP = pn->getRightMostChild();
  for ( ; child != nullptr ; child = child->getLeftSibling(), childP = childP->getLeftSibling() ) {
    backtrack(child,childP);
}
}
