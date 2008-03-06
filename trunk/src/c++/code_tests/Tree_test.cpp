//--------------------------------------------------
//                                        
// File: Tree_test.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: Tree_test.cpp,v 1.22 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------

#include "Tree.hpp"
#include <iostream>
#include <string>
#include "SequenceTree.hpp"
#include "AML_LeafLifting.hpp"
#include <math.h>
#include "SequenceTree.hpp"
#include "AML_local_improve.hpp"
#include "AML_given_edge_probabilities.hpp"

using namespace std;


static char *defaultSequences[]= {
  "a","ggaaaaaaaaaaaaaaaaaaaaaa",
  "b","gggaaaaaaaaaaaaaaaaaaaaa",
  "c","cccaaaaaaaaaaaaaaaaaaaaa",
  "d","aaagggaaaaaaaaaaaaaaaaaa",
  "e","aaaaaattaaaaaaaaaaaaaaaa",
  "f","cccgggttaaaaaaaaaaaaaaaa",
  "g","cccgggtaaaaaaaaaaaaaaaaa",
  "h","cccaagcaaaaaaaaaaaaaaaaa",
  "i","cccaagcaaaaaaaaaaaaaaaaa",
  "j","cccaagaaaaaaaaaaaaaaaaaa",
  "k","aaaaaaaaaaaaaaaaaaaaaaaa",
  "l","aaaaaaaaaaaaaaaaaaaaaaaa",
  "m","aaaaaaaaaaaaaaaaaaaaaaaa",
};  //totaly 13 pairs

typedef __gnu_cxx::hash_map<SequenceTree , int, objhash, objeq> tree2int_map;

int
main(int argc, char **argv){


  SequenceTree treeId("( (f,g)  , r, (h,i,j))");
  treeId.drawTree(cout);
  SequenceTree::NodeVector leafs;
  treeId.addLeafs(leafs);
  for(size_t i=0;i<leafs.size();i++){
    PRINT(i);
    PRINT(leafs[i]->getNodeId());
    PRINT(NAME(leafs[i]));
  }
  treeId.makeCanonical(leafs);
  treeId.drawTree(cout);


  //-------------------------------
  SEPARATOR();  SEPARATOR();  SEPARATOR();
  SequenceTree tree1("( (f,g)  , r, (h,i,j))");
  tree1.drawTree(cout);
  LINE();
  SequenceTree tree2 = tree1;
  LINE();
 
  str2int_hashmap name2id;
  tree1.createLeafNameToLeafIdMap(name2id);
  tree1.makeCanonical(name2id);
  tree2.makeCanonical(name2id);
  bool eq = (tree1 == tree2);
  LINE();PRINT(eq);
  ASSERT_EQ(eq,true);
  tree1.drawTree(cout);tree2.drawTree(cout);
  ASSERT_EQ(tree1,tree2);  
  ASSERT_EQ(tree1.equals(&tree2),true);
  ASSERT_EQ(tree2.equals(&tree1),true);
  PRINT(tree1.hashCode());PRINT(tree2.hashCode());
  ASSERT_EQ(tree1.hashCode(),tree2.hashCode());
  
  tree2 = SequenceTree ("( (f,g,h) , r, (i,j))");
  tree2.makeCanonical(name2id);
  eq = (tree1 == tree2);
  LINE();PRINT(eq);
  ASSERT_EQ(eq,false);

  tree1.createLeafNameToLeafIdMap(name2id);
  tree1.makeCanonical(name2id);
  tree2.makeCanonical(name2id);
  tree1.drawTree(cout);tree2.drawTree(cout);
  ASSERT_EQ(tree1==tree2,false);
  ASSERT_EQ(tree1.equals(&tree2),false);
  ASSERT_EQ(tree2.equals(&tree1),false);
  PRINT(tree1.hashCode());PRINT(tree2.hashCode());
  ASSERT_EQ(tree1.hashCode()==tree2.hashCode(),false);
  
  tree2 = SequenceTree ("( ((f,r),g), i,(h,j))");
  tree2.makeCanonical(name2id);
  eq = (tree1==tree2);
  LINE();PRINT(eq);
  ASSERT_EQ(eq,false);

  tree1.createLeafNameToLeafIdMap(name2id);
  tree1.makeCanonical(name2id);
  tree2.makeCanonical(name2id);
  tree1.drawTree(cout);tree2.drawTree(cout);
  ASSERT_EQ(tree1==tree2,false);
  ASSERT_EQ(tree1.equals(&tree2),false);
  ASSERT_EQ(tree2.equals(&tree1),false);
  PRINT(tree1.hashCode());PRINT(tree2.hashCode());
  ASSERT_EQ(tree1.hashCode()==tree2.hashCode(),false);
  
  tree1 = SequenceTree("((f,r),g,((h,j),i))");
  tree1.makeCanonical(name2id);
  eq = (tree1==tree2);
  LINE();PRINT(eq);
  tree1.drawTree(cout);
  tree2.drawTree(cout);
  ASSERT_EQ(eq,true);
  tree1.createLeafNameToLeafIdMap(name2id);
  tree1.makeCanonical(name2id);
  tree2.makeCanonical(name2id);
  tree1.drawTree(cout);tree2.drawTree(cout);
  ASSERT_EQ(tree1==tree2,true);
  ASSERT_EQ(tree1.equals(&tree2),true);
  ASSERT_EQ(tree2.equals(&tree1),true);
  PRINT(tree1.hashCode());PRINT(tree2.hashCode());
  ASSERT_EQ(tree1.hashCode(),tree2.hashCode());
  
  SEPARATOR();  SEPARATOR();  SEPARATOR();

  exit(1); 
  //=-=======
  SequenceTree tree("( (a:2,b:3):4  , (e:5,f:6,g:7):8)");

  cout << tree << endl;

  SequenceTree::NodeVector vec;

  
  cout << " ---------- DRAW " << endl;
  tree.drawTree(cout);

  
  cout << " ----- LEAFS " << endl;
  vec.clear();
  tree.addLeafs(vec);

  for ( size_t i = 0 ; i < vec.size() ; i++ )
    cout << vec[i] << endl;

    
  cout << " ----- POSTFIX " << endl;
  vec.clear();
  tree.addNodesInPostfixOrder(vec);

  for ( size_t i = 0 ; i < vec.size() ; i++ )
    cout << vec[i] << endl;

  
  cout << " ----- PREFIX " << endl;
  vec.clear();
  tree.addNodesInPrefixOrder(vec);

  for ( size_t i = 0 ; i < vec.size() ; i++ )
    cout << vec[i] << endl;

  

  cout << " ---------- DRAW SUBTREE" << endl;
  SequenceTree::Node *n = vec[1];
  cout << n << endl;
  tree.drawSubtree(cout,n);

  
  cout << " ----- LEAFS SUBTREE" << endl;

  vec.clear();
  tree.addLeafs(vec,n);

  for ( size_t i = 0 ; i < vec.size() ; i++ )
    cout << vec[i] << endl;

    
  cout << " ----- POSTFIX SUBTREE" << endl;
  vec.clear();
  tree.addNodesInPostfixOrder(vec,n);

  for ( size_t i = 0 ; i < vec.size() ; i++ )
    cout << vec[i] << endl;

  
  cout << " ----- PREFIX SUBTREE" << endl;
  vec.clear();
  tree.addNodesInPrefixOrder(vec,n);

  for ( size_t i = 0 ; i < vec.size() ; i++ )
    cout << vec[i] << endl;



  cout << " ---------- COPY CONSTRUCTOR" << endl;

  SequenceTree t1(tree);
  t1.drawTree(cout);

  cout << " ---------- COPY SUBTREE" << endl;
  cout << vec[1] << endl;
  SequenceTree t2(*(vec[1]));
  t2.drawTree(cout);
  
  cout << " ---------- ASSIGNMENT" << endl;
  t1 = t2;
  t2.drawTree(cout);

  cout <<" -------------- JOIN " << endl;
  t1 = tree;
  t1.joinTreeAtRoot(t2);
  t1.drawTree(cout);

  cout << "--------------- JOIN AT NODE " << endl;
  t1 = SequenceTree("( (a:2,b:3):4  , (e:5,f:6,g:7):8)");;
  t2 = SequenceTree("( (x:2,y:3):4  , (z:5,u:6,v:7):8)");;
  t1.drawTree(cout);
  t2.drawTree(cout);
  vec.clear();
  t1.addNodesInPrefixOrder(vec);
  cout <<" ***** joining tree at " << vec[1] <<endl; 
  t1.joinTreeAtNode(vec[1],t2);
  t1.drawTree(cout);

  cout << "--------------- JOIN AT NODES " << endl;
  t1 = SequenceTree("( (a:2,b:3):4  , (e:5,f:6,g:7):8)");;
  t2 = SequenceTree("( (x:2,y:3):4  , (z:5,u:6,v:7):8)");;
  t1.drawTree(cout);
  t2.drawTree(cout);
  vec.clear();
  t1.addNodesInPrefixOrder(vec);
  cout <<" ***** joining tree at " << vec[1] <<endl;
  SequenceTree::Node *firstN = vec[1];
  vec.clear();
  t2.addNodesInPrefixOrder(vec);
  cout <<" ***** and at " << vec[1] <<endl;
  SequenceTree::Node *secondN = vec[1];
  t1.joinTreeAtNodes(firstN,secondN);
  t1.drawTree(cout);
  
  cout <<"---------------- PATH " << endl;
  vec.clear();
  tree.addLeafs(vec);
  SequenceTree::Node *n1 = vec[0];
  SequenceTree::Node *n2 = vec[vec.size()-1];
  vec.clear();
  tree.drawTree(cout);
  cout << " PATH BETWEEN " << n1 << " and  " << n2 << endl;
  tree.addNodesOnPath(vec,n1,n2);
  for ( size_t i = 0 ; i < vec.size() ; i++ )
    cout << vec[i] << endl;


  cout << "--------------- REMOVE NODE " << endl;
  SequenceTree t3("( (a:2,b:3):4  , (e:5,f:6,g:7):8)");  
  t3.drawTree(cout);
  n = t3.getRoot()->getRightMostChild();
  cout << "REMOVING   " << n << endl;
  t3.removeAndDelete(n);
  t3.drawTree(cout);

  cout << "--------------- REROOT " << endl;
  
  t3 = SequenceTree("(( (a:2,b:3):4  , (e:5,f:6,g:7):8),( (x:2,y:3):4  , (z:5,h:6,i:7):8))");  
  t3.drawTree(cout);
  vec.clear();
  t3.addNodesInPostfixOrder(vec);
  cout << "*** rerooting at " << vec[0] << endl;
  t3.reRootAt(vec[0]);
  t3.drawTree(cout);

  cout << "--------------- SHORTCUT " << endl;
  t3 = SequenceTree("(((E1:5)NOTBEHERE:5,F)A, (G,H)B:3,(I,J)C:3)NOTBEHERE");
  t3.drawTree(cout);
  vec.clear();
  t3.addNodesWithDegree(vec,2);
  assert ( vec.size() == 1 );
  PRINT(vec[0]);
  t3.shortcutNode(vec[0]);
  t3.drawTree(cout);
  PRINT(t3.shortcutNode(t3.getRoot()));
  t3.drawTree(cout);

  cout << "---------------- SHORTCUT" << endl;
  t3 = SequenceTree("(((E1:5)NOTBEHERE:5,F)A:3, (I,J)C:3)NOTBEHERE");
  t3.drawTree(cout);
  PRINT_EXP(t3.shortcutDegree2Nodes());
  t3.drawTree(cout);
  
  
  cout << "FINNISHED SHORTCUTTING" << endl;
  cout << "********************************" << endl;
  cout << "********************************" << endl;
  cout << "********************************" << endl;
  cout << "Liklihood Tests " << endl;

  SequenceTree dt("( (a:0.1,b:0.4)c:0.1  , (d:0.1,e:0.2,f:0.4)g:0.2)h");
  dt.mapSequencesOntoTree(defaultSequences,dt.getNumNodes());
  dt.verbosePrint(cout);
    
  PRINT(dt.compute_loglikelihood());

  cout << "-------------------- optimal likelihood " << endl;
  vec.clear();
  dt.addNodesInPrefixOrder(vec);
  for ( size_t i = 0 ; i<vec.size() ; i++)
    vec[i]->data.dbl = -1;
  dt.verbosePrint(cout);
  PRINT(dt.compute_loglikelihood());

  
  cout << "********************************" << endl;
  cout << "********************************" << endl;
  cout << "********************************" << endl;
  cout << "LEAF LIFTING " << endl;

  dt = SequenceTree("( (a,b)  , (c,d,e))");
  dt.mapSequencesOntoTree(defaultSequences,dt.getNumNodes());
  dt.verbosePrint(cout);
  
  float likelihood;
  PRINT(likelihood = computeOptimal_AML_LeafLifting(dt, P_DISTANCE));
  dt.verbosePrint(cout);

  cout << "Compare with other likelihood " << endl;
  float tmp;
  PRINT( tmp = dt.compute_loglikelihood());
  
  assert( fabs(likelihood - tmp) <0.00002 );

  cout << "-------------------------" << endl;
  tree = SequenceTree("( (f,g)  , (h,i,j))");
  tree.mapSequencesOntoTree(defaultSequences,13);
  dt.joinTreeAtRoot(tree);
  dt.verbosePrint(cout);
  PRINT(likelihood = computeOptimal_AML_LeafLifting(dt, P_DISTANCE));
  dt.verbosePrint(cout);
  
  cout << "Compare with other likelihood " << endl;
  PRINT(tmp = dt.compute_loglikelihood());

  assert( fabs(likelihood - tmp) <0.00002 );

  PRINT_EXP(AML_local_improve(dt,P_DISTANCE));
  PRINT(dt.compute_loglikelihood());
  dt.verbosePrint(cout);
  cout << "------------------------- AML sequences given edge probabilities" << endl;
  PRINT(computeAML_given_edge_probabilities(dt));
  PRINT(dt.compute_loglikelihood());
  dt.verbosePrint(cout);
  cout << "------------------------- PARSIMONY" << endl;
  PRINT(dt.computeMostParsimoniousSequences());
  dt.verbosePrint(cout);
  PRINT(dt.compute_loglikelihood());
  PRINT_EXP(AML_local_improve(dt,P_DISTANCE));
  PRINT(dt.compute_loglikelihood());

  return 1;
}




