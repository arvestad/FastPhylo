//--------------------------------------------------
//                                        
// File: AML_local_improve.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: AML_local_improve.cpp,v 1.1 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------

#include "AML_local_improve.hpp"
#include <string>
#include "AML_star.hpp"

using namespace std;
void
AML_local_improve(SequenceTree &t, sequence_model model){

  //get all binary nodes
  SequenceTree::NodeVector vec;
  vec.reserve(t.getNumNodes());
  t.addNodesWithDegree(vec,3);


  //loop until no improvement
  string *neigh[3];
  int size = vec.size();
  for ( int i = 0 ; i < size ; i++ ){
    SequenceTree::Node *n = vec[i];

    //fetch the neighboring sequences
    SequenceTree::Node *child = n->getRightMostChild();
    neigh[0] = &(child->data.s.seq);
    child = child->getLeftSibling();
    neigh[1] = &(child->data.s.seq);
    child = child->getLeftSibling();
    if ( child != NULL )
      neigh[2] = &(child->data.s.seq);
    else
      neigh[2] = &(n->getParent()->data.s.seq);

    //compute the current likelihood
    if (improve_AML_star(*neigh[0],*neigh[1],*neigh[2], n->data.s.seq, model)){
      //PRINT(n->getNodeId());      PRINT(n->data.s.seq);
      i = -1;
      continue;
    }
  }   
}
