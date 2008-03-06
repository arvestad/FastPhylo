//--------------------------------------------------
//                                        
// File: AML_given_edge_probabilities.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: AML_given_edge_probabilities.cpp,v 1.2 2006/12/19 09:24:36 isaac Exp $                                 
//
//--------------------------------------------------

#include "AML_given_edge_probabilities.hpp"

#include <string>
#include "nucleotide.hpp"
#include <vector>
#include "log_utils.hpp"
#include <iostream>

using namespace std;

#define MAX(a,b) ((a) > (b) ? (a) : (b))

//The loglikelihood information is saved in a vector of 4 integers. Each position
//has one double representing the logliklihood of labeling with that nucleotide
typedef vector<double> p_info;

typedef std::vector<p_info> p_vector;


std::istream &
operator>>(std::istream &in,p_info p){
  in >>p[DNA_A_]; in >> p[DNA_C_];in >> p[DNA_G_];in >> p[DNA_T_];
  return in;
}

std::ostream&
operator<<(std::ostream & os, p_info  p){
  os << "a="<<p[DNA_A_] << " c="<<p[DNA_C_] << " g="<<p[DNA_G_] << " t="<<p[DNA_T_];
  return os;
}

typedef Tree<p_vector> AMLTree;

static void
backtrack(SequenceTree::Node *sn, AMLTree::Node *pn);


double
computeAML_given_edge_probabilities(SequenceTree &resultTree){

  SequenceTree::Node *oldRoot = resultTree.getRoot();
  resultTree.reRootAtHighDegreeNode();

  p_vector dummy;
  AMLTree pTree(resultTree,dummy);
  
  
  AMLTree::NodeVector nodes;
  nodes.reserve(pTree.getNumNodes());
  SequenceTree::NodeVector origNodes;
  origNodes.reserve(pTree.getNumNodes());
  
  //init the leafs
  p_info defaultVal(4);
  defaultVal[0]=-100000;defaultVal[1]=-100000;defaultVal[2]=-100000;defaultVal[3]=-100000;
  resultTree.addLeafs(origNodes);
  pTree.addLeafs(nodes);
  string &seq = oldRoot->data.s.seq;
  for ( size_t i = 0 ;i < origNodes.size() ; i++ ){
    seq = origNodes[i]->data.s.seq;
    p_vector &pos = nodes[i]->data;
   
    pos.resize(seq.length(),defaultVal);
   
    for ( size_t j = 0 ; j < seq.length() ; j++ ){
      nucleotide n = char2nucleotide(seq[j]);
      pos[j][n] = 0;
    }
    //    cout << pos << endl;
  }

  //constants
  const size_t seqlen = seq.length();
  const size_t numNodes = resultTree.getNumNodes();
  const size_t numSymbols = 4;

    
  //BOTTOM UP
  nodes.clear();
  pTree.addNodesInPostfixOrder(nodes);
  origNodes.clear();
  resultTree.addNodesInPostfixOrder(origNodes);
  
  for ( size_t i = 0 ; i < numNodes ; i++ ){
    if ( nodes[i]->isLeaf() ) continue;
    p_vector &pos = nodes[i]->data;

    //PRINT("----------------------------------");PRINT(nodes[i]->getNodeId());
    
    pos.resize(seq.length(),defaultVal);
    //for each position in the sequence
    for( size_t j = 0 ; j < seqlen ; j++ ){
      p_info &p = pos[j];
      //PRINT(j);PRINT(p);
      //For each possible symbol compute the loglikelihood score
      //by taking the maximum sum over the children 
      for ( size_t sym = 0 ; sym < numSymbols ; sym++ ){
        p[sym] = 0;
        AMLTree::Node *child = nodes[i]->getRightMostChild();
        SequenceTree::Node *s_child = origNodes[i]->getRightMostChild();
        for ( ; child != NULL ; child = child->getLeftSibling() , s_child = s_child->getLeftSibling() ){
          p_info &child_p = (child->data)[j];
          double edge_prob = s_child->data.dbl;
          double score = child_p[0] + log( sym != 0 ? edge_prob/3.0 : (1-edge_prob) );
          for ( size_t symC = 1 ; symC < numSymbols ; symC++ ){
            double symCp = child_p[symC] + log( sym != symC ? edge_prob/3.0 : (1-edge_prob) );
            score = MAX(score, symCp);
          }
          //add the best cost the loglikelihood at the parent
          //PRINT(score);
          p[sym] += score;
        }//end children loop
      }//end symbol loop
    }//end sequence position
  }//end node loop


  //BACK TRACK
  //first create the best sequence in the root and then backtrack down.
  SequenceTree::Node *sroot = resultTree.getRoot();
  AMLTree::Node *proot = pTree.getRoot();
  p_vector &pos = proot->data;
  seq = sroot->data.s.seq;
  seq.resize(seqlen,'x');
  double totalLoglikelihood = 0;
  for ( size_t i = 0 ; i < seqlen ; i++ ){
    p_info &p = pos[i];
    nucleotide max_n = DNA_A_;
    double score = p[DNA_A_];
    if ( p[DNA_C_] > score ){
      max_n = DNA_C_;score = p[DNA_C_];
    }
    if ( p[DNA_G_] > score ){
      max_n = DNA_G_;score = p[DNA_G_];
    }
    if ( p[DNA_T_] > score ){
      max_n = DNA_T_;score = p[DNA_T_];
    }
    seq[i] = nucleotide2char(max_n);
    totalLoglikelihood += score;
  }
  SequenceTree::Node *child = sroot->getRightMostChild();
  AMLTree::Node *childP = proot->getRightMostChild();
  for ( ; child != NULL ; child = child->getLeftSibling(), childP = childP->getLeftSibling() )
    backtrack(child,childP);
  
  //reroot the tree to its old position
  resultTree.reRootAt(oldRoot);

  return totalLoglikelihood;
}



static void
backtrack(SequenceTree::Node *sn, AMLTree::Node *pn){

  if ( sn->isLeaf() )
    return;

  
  p_vector &pos = pn->data;
  const size_t seqlen = pos.size();
  string &seq = sn->data.s.seq;
  seq.resize(seqlen,'x');
  string &parentseq = sn->getParent()->data.s.seq;
  double edge_prob = sn->data.dbl;
  
  for ( size_t i = 0 ; i < seqlen ; i++ ){
    p_info &p = pos[i];
    nucleotide parent_n = char2nucleotide(parentseq[i]);
    nucleotide max_n = DNA_A_;
    double score = p[DNA_A_] + (parent_n == DNA_A_ ? log(1-edge_prob) : log(edge_prob/3.0));
    double tmpscore = p[DNA_C_] + (parent_n == DNA_C_ ? log(1-edge_prob) : log(edge_prob/3.0));
    if ( tmpscore > score ){
      max_n = DNA_C_;
      score = tmpscore;
    }
    tmpscore = p[DNA_G_] + (parent_n == DNA_G_ ? log(1-edge_prob) : log(edge_prob/3.0));
    if ( tmpscore > score ){
      max_n = DNA_G_;
      score = tmpscore;
    }
    tmpscore = p[DNA_T_] + (parent_n == DNA_T_ ? log(1-edge_prob) : log(edge_prob/30.0));
    if ( tmpscore > score ){
      max_n = DNA_T_;
    }
    seq[i] = nucleotide2char(max_n);
  }

  

  //continue backtracking
  SequenceTree::Node *child = sn->getRightMostChild();
  AMLTree::Node *childP = pn->getRightMostChild();
  for ( ; child != NULL ; child = child->getLeftSibling(), childP = childP->getLeftSibling() )
    backtrack(child,childP);
}


//----------------------------
void
convert_edge_lengths_to_probabilities(SequenceTree &t, sequence_model model){

  if ( model == P_DISTANCE )
    return ;

  if ( model != JC )
    PROG_ERROR("unimplemented model");


  SequenceTree::NodeVector nodes;
  nodes.reserve(t.getNumNodes());
  t.addNodesInPrefixOrder(nodes);

  //skip the root
  for ( size_t i = 1 ; i < nodes.size() ; i++ ){

    SequenceTree::Node *n = nodes[i];
    double dist = n->data.dbl;
    //invert the JC formula
    double prob = 3.0/4.0 * ( 1 - exp(-4.0/3.0 * dist ));
    n->data.dbl = prob;
  }
}
