//--------------------------------------------------
//                                        
// File: AML_SpanningTree_test.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: Big_AML_test.cpp,v 1.3 2006/12/05 09:08:29 isaac Exp $                                 
//
//--------------------------------------------------

#include "Big_AML.hpp"
#include "AML_LeafLifting.hpp"


#include <vector>
#include <string>
#include <iostream>

#include "log_utils.hpp"
#include "dna_pairwise_sequence_likelihood.hpp"

using namespace std;




int
main ( int argc, char **argv ){


  vector<Sequence> seqs;
  seqs.push_back(Sequence("a","aaaaaaaaaaaaaaaaaaaaaaaa"));
  seqs.push_back(Sequence("b","ggggaaaaaaaaaaaaaaaaaaaa"));
  seqs.push_back(Sequence("c","ggggggggaaaaaaaaaaaaaaaa"));
  seqs.push_back(Sequence("d","ggggggggttttaaaaaaaaaaaa"));
  seqs.push_back(Sequence("e","gggggaaaaaaaaggaaaaaaaaa"));
  seqs.push_back(Sequence("f","gggggaaaaaaaaggttaaaaaaa"));
  seqs.push_back(Sequence("g","gggggaaaaaaaaggaaaaggaaa"));

  SequenceTree tree;

  PRINT_EXP(Big_AML(seqs,tree, P_DISTANCE));
  tree.verbosePrint(cout);
  PRINT(tree.compute_loglikelihood());
  tree.drawTree(cout);
  cout << "--------------------------  JC" << endl;
  
  PRINT_EXP(Big_AML(seqs,tree, JC));
  tree.verbosePrint(cout);
  PRINT(tree.compute_loglikelihood());
  tree.drawTree(cout);
  return 0;
}



