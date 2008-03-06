//--------------------------------------------------
//                                        
// File: AML_SpanningTree.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: AML_SpanningTree.cpp,v 1.1 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------

#include "AML_SpanningTree.hpp"
#include <string>
#include <stdlib.h>

#include "LikelihoodMatrix.hpp"
#include <iostream>
#include "dna_pairwise_sequence_likelihood.hpp"
#include "string_compare.hpp"

#include "NeighborJoining.hpp"

using namespace std;

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
AML_SpanningTree(std::vector<std::string *> &sequences,
                 SequenceTree &resultTree){

  //1. compute all pairwise likelihoods
  size_t numseqs = sequences.size();
  LikelihoodMatrix lm(numseqs);
  fillLikelihoodMatrix(sequences,lm);
  //  cout << lm << endl;
  
  //2. Sort all pairwise relations
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

  //3. Create a tree for each sequence
  vector<SequenceTree> trees(numseqs);
  vector<SequenceTree::Node *> nodes(numseqs);
  for ( size_t i = 0 ; i < numseqs ; i++ ){
    string_double strflt;
    strflt.str.append(*sequences[i]);
    strflt.flt = -1;
    trees[i] = SequenceTree(strflt);
    nodes[i] = trees[i].getRoot();
  }

  //4. Build the MST
  double likelihood = 0;
  size_t numjoinsperformed = 0;
  for ( size_t i = numvals-1 ; numjoinsperformed < numseqs-1 ; i-- ){
    pairwise_relation r = vals[i];
    //check that the nodes aren't in the same tree
    if ( nodes[r.i]->getTree() == nodes[r.j]->getTree() )
      continue;

    //else join the trees at the respective nodes.
    nodes[r.i]->getTree()->joinTreeAtNodes(nodes[r.i],nodes[r.j]);
    numjoinsperformed++;
    likelihood += r.likelihood;
  }

  
  //5. Set the result tree
  resultTree = *(nodes[0]->getTree());  
  return likelihood;
}



void
addLeafChildToInternalNodes(SequenceTree &resultTree){
  SequenceTree::NodeVector vec;

  resultTree.addNodesInPrefixOrder(vec);

  for ( size_t i = 0 ; i < vec.size() ; i++ ){

    SequenceTree::Node *n = vec[i];
    if ( n->getDegree() > 1 ){
      n->addChild(n->data);
    }
  }
}


typedef DistanceMatrix<SequenceTree::Node *,double,Data_init<SequenceTree::Node *>,Data_printOn<SequenceTree::Node *>,Data_init<double>,Data_printOn<double> > NJMatrix;

void
resolveHighDegreeInternalNodes(SequenceTree &resultTree){

  SequenceTree::NodeVector allnodes;
  SequenceTree::NodeVector neighboringNodes;
  NJMatrix njm(resultTree.getNumNodes());
  SequenceTree::Node *oldRoot = resultTree.getRoot();
  
  resultTree.addNodesInPrefixOrder(allnodes);
  for ( size_t node_i = 0 ; node_i < allnodes.size() ; node_i++ ){
    SequenceTree::Node *n = allnodes[node_i];
    if ( n->getDegree() > 3 ){
      resultTree.reRootAt(n);
      
      neighboringNodes.clear();
      SequenceTree::Node *child = n->getRightMostChild();
      while (child != NULL){
        //PRINT(child->getParent());
        neighboringNodes.push_back(child);
        child = child->getLeftSibling();
      }

      //create distance matrix for neighboring nodes
      njm.resize(neighboringNodes.size());
      for ( size_t i = 0 ; i < neighboringNodes.size() ; i++ ){
        njm.setIdentifier(i,neighboringNodes[i]);
        for ( size_t j = i+1 ; j < neighboringNodes.size() ; j++ ){
          TN_string_distance tndist = TN_string_compare(neighboringNodes[i]->data.str, neighboringNodes[j]->data.str);
          simple_string_distance sd = convert_TN_string_distance_to_simple(tndist);
          njm.setDistance(i,j,compute_JC(neighboringNodes[i]->data.str.size(),sd).distance);
        }
      }

      //do neighbour joining of the star and set the data of n to the
      //newly created nodes
      computeNeighborJoiningTree(njm, n->data);
    }
  }

  //PRINT(oldRoot);
  resultTree.reRootAt(oldRoot);
} 


