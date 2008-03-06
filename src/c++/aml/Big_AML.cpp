//--------------------------------------------------
//                                        
// File: Big_AML.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: Big_AML.cpp,v 1.1 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------

#include "Big_AML.hpp"
#include <string>
#include <stdlib.h>
#include "LikelihoodMatrix.hpp"
#include "DistanceMethodMatrix.hpp"
#include "NeighborJoining.hpp"
#include <iostream>
#include "AML_local_improve.hpp"

using namespace std;

void
Big_AML(std::vector<Sequence> &seqs, SequenceTree &resultTree, sequence_model model){

  //Compute the likelihoods
  size_t numseqs = seqs.size();
  SequenceTree::NodeMatrix lm(numseqs);
  std::vector<std::string *> strvec(numseqs);
  for ( size_t i = 0 ; i < numseqs ; i++ ){
    strvec[i] = &(seqs[i].seq);
  }
  fillLikelihoodMatrix(strvec, lm,model);
  
  //set the identifiers in the liklihood matrix to be unique trees
  //containing only one node
  vector<SequenceTree> trees;
  trees.reserve(numseqs);
  for ( size_t i = 0 ; i < numseqs ; i++ ){
      Sequence_double strdbl;
      strdbl.s = seqs[i];
      strdbl.dbl = -1;
      trees.push_back(SequenceTree(strdbl));
      lm.setIdentifier(i, trees[i].getRoot());
  }

  

  //Create spaning tree
  AML_SpanningTree(lm);

  //add leafs to internal nodes and clear names of the internal node
  resultTree = *((SequenceTree*)lm.getIdentifier(0)->getTree());
  SequenceTree::NodeVector nodes;
  nodes.reserve(2*numseqs);
  resultTree.addNodesInPrefixOrder(nodes);
  addLeafChildToInternalNodesInVector(nodes);

  //-------
  // NJ RESOLVE OF HIGH DEGREE NODES
  nodes.clear();
  resultTree.addNodesInPrefixOrder(nodes);
  SequenceTree::NodeVector starnodes;
  starnodes.reserve(2*numseqs);

  //for each high degree node add its neighbors to the starnodes vector
  //and call NJ
  for ( size_t node_i = 0 ; node_i < nodes.size() ; node_i++ ){
    SequenceTree::Node *n = nodes[node_i];
    if ( n->getDegree() > 3 ){
      resultTree.reRootAt(n);

      //      resultTree.drawTree(cout);
      starnodes.clear();
      SequenceTree::Node *child = n->getRightMostChild();
      while (child != NULL){
        starnodes.push_back(child);
        child = child->getLeftSibling();
      }

      //create distance matrix for neighboring nodes
      lm.resize(starnodes.size());
      strvec.resize(starnodes.size());
      for ( size_t i = 0 ; i < starnodes.size() ; i++ ){
        lm.setIdentifier(i,starnodes[i]);
        strvec[i] = &(starnodes[i]->data.s.seq); 
      }
      fillDistanceMethodMatrix(strvec,lm,model);
      fixFactorCorrection_ToMAX(lm);
      //do neighbour joining of the star and set the data of n to the
      //newly created nodes
      computeNeighborJoiningTree(lm, n->data);
    }    
  }
  //END NJ RESOLVE
  //----
  // HEURISTIC IMPROVE
  AML_local_improve(resultTree,model);
}



//-------------------------------------------------
// SPANING TREE

//contains the 
struct pairwise_relation{
  size_t i;
  size_t j;
  double likelihood;
};

static int
compare_relations(const void *v1, const void *v2){
  pairwise_relation *r1 = (pairwise_relation *)v1;
  pairwise_relation *r2 = (pairwise_relation *)v2;

  if ( r1->likelihood < r2->likelihood )
    return -1;
  if ( r1->likelihood > r2->likelihood )
    return 1;
  else
    return 0;
}

double
AML_SpanningTree(SequenceTree::NodeMatrix &lm){

  //1. Sort all pairwise relations
  size_t numseqs = lm.getSize();
  const size_t numvals = numseqs*(numseqs-1)/2;
  pairwise_relation vals[numvals];
  size_t rel_i = 0;
  
  for ( size_t i = 0 ; i < numseqs ; i++ ){
    for ( size_t j = i+1 ; j < numseqs ; j++ ){
      vals[rel_i].i = i;
      vals[rel_i].j = j;
      vals[rel_i].likelihood = lm.getDistance(i,j);
      rel_i++;
    }
  }
  assert ( rel_i == numvals );
  qsort(vals, numvals, sizeof(pairwise_relation),  &compare_relations);


  //2. Build the MST
  double likelihood = 0;
  size_t numjoinsperformed = 0;
  for ( int i = numvals-1 ; numjoinsperformed < numseqs-1 && i != -1 ; i-- ){
    pairwise_relation r = vals[i];
    //check that the nodes aren't in the same tree
    if ( lm.getIdentifier(r.i)->getTree() == lm.getIdentifier(r.j)->getTree() )
      continue;

    //else join the trees at the respective nodes.
    lm.getIdentifier(r.i)->getTree()->joinTreeAtNodes(lm.getIdentifier(r.i),lm.getIdentifier(r.j));
    numjoinsperformed++;
    likelihood += r.likelihood;
  }

  
  //3. return the sum of log likelihoods
  return likelihood;
}

//-----------------------------------------------

void
addLeafChildToInternalNodesInVector(std::vector<SequenceTree::Node *> &nodes){
  
  for ( size_t i = 0 ; i < nodes.size() ; i++ ){

    SequenceTree::Node *n = nodes[i];
    if ( n->getDegree() > 1 ){
      n->addChild(n->data);
      n->data.s.name.clear();
    }
  }
}

