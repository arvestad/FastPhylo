//--------------------------------------------------
//                                        
// File: AML_LeafLifting.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: AML_LeafLifting.cpp,v 1.2 2006/12/19 09:24:36 isaac Exp $                                 
//
//--------------------------------------------------

#include "AML_LeafLifting.hpp"
#include <iostream>
#include "stl_utils.hpp"
#include <string>
#include <float.h>
#include "LikelihoodMatrix.hpp"

using namespace std;

typedef std::vector<double> double_vector;

//template std::ostream& operator<< <double> (std::ostream &, std::vector<double> &);

typedef Tree<double_vector> LiftingTree;



static void backtrack(SequenceTree::Node *sn, LiftingTree::Node *n, size_t parentMax);

SequenceTree::NodeVector origleafs;
double MINVAL;
SequenceTree::NodeMatrix lm(10);

double
computeOptimal_AML_LeafLifting(SequenceTree &t, sequence_model model){

  //make sure that the tree is rooted at a highdegree node
  //at the end the tree is returned to its old form
  SequenceTree::Node *oldRoot = t.getRoot();
  t.reRootAtHighDegreeNode();

  const size_t numLeafs = t.getNumLeafs();
  origleafs.reserve(numLeafs);
  origleafs.clear();
  t.addLeafs(origleafs);

  
  //------------------------
  //precompute all pairwise likelihodds
  lm.resize(numLeafs);
  vector<string *> strvec(numLeafs);
  for ( size_t i = 0 ; i < numLeafs ; i++ )
    strvec[i] = &(origleafs[i]->data.s.seq);

  fillLikelihoodMatrix(strvec,lm, model);
  //  cout << lm << endl;

  //find a very small value to init the leafs with
  MINVAL = 0.0;
  for ( size_t i = 0 ; i < numLeafs ; i++ ){
    for ( size_t j = i+1 ; j < numLeafs ; j++ ){
      double tmp = lm.getDistance(i,j);
      if ( tmp < MINVAL ) MINVAL = tmp;
    }
  }
  MINVAL -=1;
  MINVAL *=5;
  
  
  //-------------------------
  // LIFTING TREE
  // the tree in which we perform the bottom up algo
  double_vector dummy;
  LiftingTree lt(t,dummy);
  
  //likelihoods[i] is the optimal likelihood of labeling the node
  //with the same labels as origleafs[i]
  
  //------------------------
  //INIT
  //init the leafs
  LiftingTree::NodeVector liftingleafs;
  lt.addLeafs(liftingleafs);
  for ( size_t i = 0 ; i < numLeafs ; i++ ){
    LiftingTree::Node *n = liftingleafs[i];
    double_vector &likelihoods = n->data;
    likelihoods.resize(numLeafs);

    for ( size_t j = 0 ; j < numLeafs ; j++ )
      likelihoods[j] = MINVAL;

    likelihoods[i] = 0.0;//only the original leaf label is allowed
  }

  //------------------------
  //BOTTOM UP
  //Do the bottom up of the internal nodes
  LiftingTree::NodeVector lnodes;
  lt.addNodesInPostfixOrder(lnodes);
  for ( size_t i = 0 ; i < lnodes.size() ; i++ ){
    LiftingTree::Node *n = lnodes[i];    
    double_vector &likelihoods = n->data;
    likelihoods.resize(numLeafs);

    if ( n->isLeaf() )
      continue;

    // likelihoods[j] <- ForAll children SUM max{ L(childlabel) + L(string_j,childlabel) }
    for ( size_t j = 0 ; j < numLeafs ; j++ ){
      likelihoods[j] = 0;
      
      LiftingTree::Node *child = n->getRightMostChild();
      for ( ; child != NULL ; child = child->getLeftSibling() ){
        double maxVal = child->data[0] + lm.getDistance(j,0);
        for ( size_t k = 1 ; k < numLeafs ; k++ ){
          double tmp = child->data[k] + lm.getDistance(j,k);
          maxVal = (maxVal > tmp ? maxVal : tmp);
        }
        likelihoods[j] += maxVal;
      }
    }    
  }

  //------------------------
  // BACKTRACK
  // Do the backtracking and set the correct leaf labels in
  // the input tree.
  SequenceTree::NodeVector orignodes;
  t.addNodesInPrefixOrder(orignodes);
  lnodes.clear();
  lt.addNodesInPrefixOrder(lnodes);

  //find max in root
  double_vector &likelihoods = lnodes[0]->data;
  size_t maxi = 0;
  double maxlikelihood = likelihoods[0];
  for ( size_t i = 1 ; i < numLeafs ; i++ ){
    if ( likelihoods[i] > maxlikelihood ){
      maxi = i;
      maxlikelihood = likelihoods[i];
    }
  }

  orignodes[0]->data.s.seq = origleafs[maxi]->data.s.seq;

  //backtrack down
  SequenceTree::Node *sn_child = orignodes[0]->getRightMostChild();
  LiftingTree::Node *n_child = lnodes[0]->getRightMostChild();
  while ( n_child != NULL ){
    backtrack(sn_child,n_child,maxi);
    n_child = n_child->getLeftSibling();
    sn_child = sn_child->getLeftSibling();
  }
  

  //-----------------------
  //reroot at old root
  t.reRootAt(oldRoot);
  return maxlikelihood;
}



static void
backtrack(SequenceTree::Node *sn, LiftingTree::Node *n, size_t parentMax){

  if ( sn->isLeaf() ){
    return;
  }

  size_t maxi = 0;
  double maxVal = n->data[0]+ lm.getDistance(parentMax,0);
  for ( size_t i = 0 ; i < origleafs.size() ; i++ ){
    double tmp = n->data[i]+ lm.getDistance(parentMax,i);
    if ( tmp > maxVal ){
      maxi = i;
      maxVal = tmp;
    }
  }

  sn->data.s.seq = origleafs[maxi]->data.s.seq;
  sn->data.dbl = -1;

  SequenceTree::Node *sn_child = sn->getRightMostChild();
  LiftingTree::Node *n_child = n->getRightMostChild();
  while ( n_child != NULL ){
    backtrack(sn_child,n_child,maxi);
    n_child = n_child->getLeftSibling();
    sn_child = sn_child->getLeftSibling();
  }
}
