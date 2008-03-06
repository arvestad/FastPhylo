//--------------------------------------------------
//                                        
// File: SequenceBasedNJ.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: SequenceBasedNJ.cpp,v 1.2 2006/04/23 08:06:50 isaac Exp $                                 
//
//--------------------------------------------------

#include "SequenceBasedNJ.hpp"
#include <string>
#include "DNA_b128_String.hpp"
#include "log_utils.hpp"
#include "Tree.hpp"
#include <iostream>
#include <float.h>


//A tree of DNA_b128_String s

struct empty_b128_init{void operator()(std::istream &in, DNA_b128_String &s) const{}};
typedef Tree<DNA_b128_String,empty_b128_init, Data_printOn<DNA_b128_String> > b128Tree;

//A DistanceMatrix with identifiers of b128Tree::Nodes
struct empty_b128Tree_init{void operator()(std::istream &in, b128Tree::Node &s) const{}};

typedef DistanceMatrix< b128Tree::Node *, double, empty_b128Tree_init, Data_printOn<b128Tree::Node *>,
                        Data_init<double>, Data_printOn<double> > b128Matrix;

//----------------------------
//FORWARD DECLARTIONS
double
computeK2PDistance(DNA_b128_String &s1, DNA_b128_String &s2);
void
fillMatrix(b128Matrix &dm);


//--------------------------------------------------
// THE SEQUENCE BASED NJ ALGO
//
//
void
computeSequenceBasedNJ(std::vector<Sequence> &seqs, SequenceTree &resultTree){

  // 1. Create a star tree with the leafs being the input sequences in b128 format.
  DNA_b128_String defaultString(seqs[0].seq.size());
  b128Tree tree(defaultString);
  b128Tree::Node *root = tree.getRoot();
  obj_ptr2obj_ptr_hashmap node2seqs((size_t)(seqs.size()*1.5));
  b128Tree::NodeVector leafs;
  for ( size_t i = 0 ; i < seqs.size()  ; i++ ){
    b128Tree::Node *leaf = root->addChild(defaultString);
    node2seqs[leaf] = &(seqs[i]);
    leafs.push_back(leaf);
    (leaf->data).append(seqs[i].seq);
  }
  
  // 2. Compute the DistanceMatrix for the seqs.
  b128Matrix dm(seqs.size());
  for ( size_t i = 0 ; i < leafs.size() ; i++ ){
    dm.setIdentifier(i,leafs[i]);
  }

  fillMatrix(dm);

  //std::cout << dm << std::endl;
  // 3. COMPUTE ROW SUMS
  double rowSums[dm.getSize()];
  for ( size_t row = 0 ; row < dm.getSize() ; row++ ){
    double sum = 0;
    size_t i =0;
    for ( ; i < dm.getSize() ; i++ )
      sum += dm.getDistance(row,i);

    rowSums[row] = sum;
  }

  //----------------
  // 4. 
  // NJ ITERATION
  //compute the row sums

  while ( dm.getSize() > 3 ) {

    //FIND MIN PAIR
    //find the minimal value
    double minVal = FLT_MAX;
    size_t mini = 1000000;
    size_t minj = 1000000;
    for ( size_t i = 0 ; i < dm.getSize() ; i++ ){
      for ( size_t j = i+1 ; j < dm.getSize() ; j++ ){
        double newVal = (dm.getSize() - 2.0)*dm.getDistance(i,j) - rowSums[i] - rowSums[j];
        //std::cout << newVal << " , ";
        if ( newVal < minVal ){
          minVal = newVal;
          mini = i;
          minj = j;
        }
      }
    }
    //    std::cout << std::endl;
    //PRINT(minVal);
    //make sure that minj is the last row in the matrix
    if ( mini == dm.getSize() -1 ){
      mini = minj;
    }
    else {
      dm.swapRowToLast(minj);
      double tmp = rowSums[dm.getSize()-1];
      rowSums[dm.getSize()-1] = rowSums[minj];
      rowSums[minj] = tmp;
    }
    minj = dm.getSize()-1;

    //CLUSTER THE LEAFS
    DNA_b128_String &child1str = dm.getIdentifier(mini)->data;
    DNA_b128_String &child2str = dm.getIdentifier(minj)->data;
    b128Tree::Node *parent = dm.getIdentifier(mini)->getTree()->detachFromParentAndAddAsSiblings(dm.getIdentifier(mini),dm.getIdentifier(minj), defaultString);
    dm.setIdentifier(mini, parent);

    //COMPUTE PARSIMONY AND SET IN PARENT
    DNA_b128_String &parentstr = parent->data;
    DNA_b128_String::create_weighted_parsimonious(parentstr,child1str,child2str);
    //COMPUTE DISTANCES FROM PARENT TO ALL OTHER NODES
    //PRINT(mini);PRINT(minj);
    for ( size_t i = 0 ; i < dm.getSize()-1 ; i++ ){//skip last row
      double dist2iandj = dm.getDistance(mini,i) + dm.getDistance(minj,i);
      DNA_b128_String &leafstr = dm.getIdentifier(i)->data;
      double dist = computeK2PDistance(parentstr,leafstr);
      // regular nj update function:
      //double regnj = dist2iandj * 0.5; 
      //double studier = (dist2iandj-dm.getDistance(mini,minj))*0.5;
      //PRINT(dist); PRINT(regnj);PRINT(dist2iandj);PRINT(dist-regnj);PRINT(dist-studier);
      //PRINT(dist - dm.getDistance(mini,i) );PRINT(dist - dm.getDistance(minj,i) );
      
      dm.setDistance(mini,i, dist); 
    
      //update rowsums
      rowSums[i] = rowSums[i] - dist2iandj + dm.getDistance(mini,i);
      //PRINT(rowSums[i]);
    }

    dm.setDistance(mini,mini,0);
    
    //remove the last row of the matrix
    dm.removeLastRow();
    
    //recompute the row sum for the parent
    double sum = 0;
    for ( size_t i = 0 ; i < dm.getSize() ; i++ )
      sum += dm.getDistance(mini,i);
    rowSums[mini] = sum;

  }
  //END ITERATION
  //----------------------------------

  
  //CONVERT THE TREE TO A SEQUENCE TREE
  tree.recalcNodeStructure();
  //  tree.drawTree(std::cout);
  b128Tree::NodeVector leafnodes;
  tree.addLeafs(leafnodes);

  Sequence_double dummy;
  dummy.dbl = -1;
  resultTree = SequenceTree(tree,dummy);

  SequenceTree::NodeVector seqnodes;
  resultTree.addLeafs(seqnodes);
  
  for ( size_t i = 0 ; i < seqnodes.size() ; i++ ){
    seqnodes[i]->data.s = *((Sequence *) node2seqs[leafnodes[i]]);
  }
  //resultTree.drawTree(std::cout);
}


//
//
double
computeK2PDistance(DNA_b128_String &s1, DNA_b128_String &s2){
   simple_string_distance sd = DNA_b128_String::computeDistance(s1,s2);
   ML_string_distance ml_dist = compute_K2P_fixratio(s1.getNumChars(),sd,2.0);
   sd = DNA_b128_String::correctDistanceWithAmbiguitiesUsingTransitionProbabilities(sd,ml_dist,s1,s2);
   ml_dist = compute_K2P_fixratio(s1.getNumChars(),sd,2.0);

   return ml_dist.distance;

   //THIS is a little bit less accurate
   //    simple_string_distance sd = DNA_b128_String::computeDistance(s1,s2);
   //    sd = DNA_b128_String::correctDistanceWithAmbiguities(sd,s1,s2);
   //    return compute_K2P_fixratio(s1.getNumChars(),sd,2.0).distance;
}

void
fillMatrix(b128Matrix &dm){
  int  numSequences = dm.getSize();
  
  for ( int i = 0 ; i < numSequences ; i++ ){
    dm.setDistance(i,i,0);
    DNA_b128_String &si = dm.getIdentifier(i)->data;
    for ( int j = i+1 ; j < numSequences ; j++ ){
      dm.setDistance(i,j,computeK2PDistance(si,dm.getIdentifier(j)->data));
    }
  }

}
