
//--------------------------------------------------
//
// File: SequenceTree.hpp
//
// Author: Mehmood Alam Khan, Isaac Elias
// e-mail: malagori@kth.se, isaac@nada.kth.se
//
//
//
//--------------------------------------------------
#include <math.h>
#include "SequenceTree.hpp"
#include "LeastSquaresFit.hpp"

using namespace std;

#define ID(node) (node)->getNodeId()

double
computeLeastSquaresEdgeLengths(const StrDblMatrix &orig_dm,  SequenceTree &tree){
  

  StrDblMatrix dm(orig_dm);
  const int numOriginalLeafs = dm.getSize();
  SequenceTree::NodeVector nodes;
  tree.recalcNodeIdsPostfixOrderAndAddInOrder(nodes);
  size_t nodeIdToRowIndex[nodes.size()];
  size_t rowIndexToNodeId[nodes.size()];
  str2int_hashmap name2Id((int)(nodes.size()*1.7));

  for(size_t i=0 ; i<nodes.size() ; i++)
    if(nodes[i]->isLeaf()){
      //PRINT(NAME(nodes[i]));PRINT(ID(nodes[i]));
      name2Id[NAME(nodes[i])] = ID(nodes[i]);
    }

  for(size_t row=0 ; row<dm.getSize() ; row++){
    str2int_hashmap::iterator f = name2Id.find(dm.getIdentifier(row));
    if(f==name2Id.end())
      USER_ERROR("name doesn't exist in tree: " << dm.getIdentifier(row));
    
    nodeIdToRowIndex[(*f).second] = row;
    rowIndexToNodeId[row] = (*f).second;
  }
  
  //the number of leafs below each node
  int numNodesBelow[nodes.size()];
  for(size_t i=0;i<nodes.size();i++)
    numNodesBelow[i]=1;

  //--------------------------------
  //BOTTOM UP TRAVERSAL IN TREE
  for(size_t i=0;i<nodes.size()-1;i++){
    if(nodes[i]->isLeaf())
      continue;

    //get the children and do the UNJ calculation to get the edge lengths
    SequenceTree::Node *parent = nodes[i];
    SequenceTree::Node *child1 = parent->getRightMostChild();
    SequenceTree::Node *child2 = child1->getLeftSibling();
    if(child2->getLeftSibling()!=NULL ){
      USER_ERROR("Have to be unrooted binary tree. Parent has " << parent->getNumChildren() << " children");
    }
    numNodesBelow[ID(parent)] = numNodesBelow[ID(child1)] + numNodesBelow[ID(child2)];
    //SEPARATOR();PRINT(NAME(child1));PRINT(NAME(child2));

    
    double sum = 0;
    for(size_t row=0;row<dm.getSize();row++){
      if(row==nodeIdToRowIndex[ID(child1)] || 
	 row==nodeIdToRowIndex[ID(child2)] )
	continue;
      sum += numNodesBelow[rowIndexToNodeId[row]]*(dm.getDistance(nodeIdToRowIndex[ID(child1)],row)-
						   dm.getDistance(nodeIdToRowIndex[ID(child2)],row));				
    }

    if(!isfinite(sum)){
      USER_ERROR("Distance Matrix contains a non finite number: " << sum);
    }

    EDGE(child1) = 0.5*dm.getDistance(nodeIdToRowIndex[ID(child1)],
				      nodeIdToRowIndex[ID(child2)])
      + 1.0/(2*(numOriginalLeafs-numNodesBelow[ID(parent)]))*sum;
    EDGE(child2) = 0.5*dm.getDistance(nodeIdToRowIndex[ID(child1)],
				      nodeIdToRowIndex[ID(child2)])
      - 1.0/(2*(numOriginalLeafs-numNodesBelow[ID(parent)]))*sum;
    // PRINT(dm.getDistance(nodeIdToRowIndex[ID(child1)],nodeIdToRowIndex[ID(child2)]));
    //     PRINT((numOriginalLeafs-numNodesBelow[ID(parent)]));
    //     PRINT(sum);PRINT(EDGE(child1));PRINT( EDGE(child2));
    //PRINT(1/(2*(numOriginalLeafs-numNodesBelow[ID(parent)]))*sum);
    
    //swap child1 to last row 
    int idOnLastRow = rowIndexToNodeId[dm.getSize()-1];
    if(idOnLastRow!=ID(child1)){
      int rowChild1 = nodeIdToRowIndex[ID(child1)];
      //PRINT(nodeIdToRowIndex[ID(child1)]);PRINT(dm.getSize());
      dm.swapRowToLast(nodeIdToRowIndex[ID(child1)]);
      nodeIdToRowIndex[idOnLastRow] = rowChild1;
      rowIndexToNodeId[rowChild1] = idOnLastRow;
      rowIndexToNodeId[dm.getSize()-1] = ID(child1);
      nodeIdToRowIndex[ID(child1)] = dm.getSize()-1;
    }
    //update distances to parent
    double w1 = (1.0*numNodesBelow[ID(child1)])/numNodesBelow[ID(parent)];
    double w2 = (1.0*numNodesBelow[ID(child2)])/numNodesBelow[ID(parent)];
    double distChild1Child2 = w1*EDGE(child1)+w2*EDGE(child2);
  
    //put parent on the row of child 2
    nodeIdToRowIndex[ID(parent)] = nodeIdToRowIndex[ID(child2)];
    rowIndexToNodeId[nodeIdToRowIndex[ID(parent)]] = ID(parent);
    int parentRow = nodeIdToRowIndex[ID(parent)]; 
    int child1Row = nodeIdToRowIndex[ID(child1)];
    int child2Row = nodeIdToRowIndex[ID(child2)];

    for(size_t row=0 ; row<dm.getSize()-1 ; row++){
      dm.setDistance(parentRow,row,
		     w1*dm.getDistance(child1Row,row)+
		     w2*dm.getDistance(child2Row,row)-
		     distChild1Child2);
    }    
    
    dm.setDistance(nodeIdToRowIndex[ID(parent)],nodeIdToRowIndex[ID(parent)],0.0);
    //remove last row
    dm.removeLastRow();
  }

  //Take care of root
  SequenceTree::Node *root = nodes[nodes.size()-1];
  if(!root->isRoot() || root->getNumChildren()!=3){
    USER_ERROR("Have to be unrooted binary tree. Root has " << root->getNumChildren() << " children");
  }
  
  //  cout << dm << endl;
  SequenceTree::Node *c1 = root->getRightMostChild();
  SequenceTree::Node *c2 = c1->getLeftSibling();
  SequenceTree::Node *c3 = c2->getLeftSibling();
  
  int c1row = nodeIdToRowIndex[ID(c1)];
  int c2row = nodeIdToRowIndex[ID(c2)];
  int c3row = nodeIdToRowIndex[ID(c3)];

  EDGE(c1) = 0.5*(dm.getDistance(c1row,c2row) + dm.getDistance(c1row,c3row)-dm.getDistance(c2row,c3row));
  EDGE(c2) = 0.5*(dm.getDistance(c2row,c1row) + dm.getDistance(c2row,c3row)-dm.getDistance(c1row,c3row));
  EDGE(c3) = 0.5*(dm.getDistance(c3row,c2row) + dm.getDistance(c3row,c1row)-dm.getDistance(c2row,c1row));
  EDGE(root) = 0;  

  //COMPUTE THE L2 SCORE
  StrDblMatrix treeM(tree.getNumLeafs());
  tree.tree2distanceMatrix(treeM);
  return computeL2(treeM, orig_dm);
}


//--------------------------------------

double
computeL2(const StrDblMatrix &A,  const StrDblMatrix &B){
  
  ASSERT_EQ(A.getSize(),B.getSize());
  
  str2int_hashmap name2Row((int)(A.getSize()*1.7));
  for(size_t i=0 ; i<A.getSize() ; i++){
    name2Row[A.getIdentifier(i)] = i;
  }

  double l2sum = 0;
  for(size_t i=0 ; i<B.getSize() ; i++){
    for(size_t j=i+1 ; j<B.getSize() ; j++ ){
      double b = B.getDistance(i,j);
      double a = A.getDistance(name2Row[B.getIdentifier(i)],
			       name2Row[B.getIdentifier(j)]);
      l2sum += (a-b)*(a-b);
    }
  }

  return l2sum;
}

// ------------ warka dang----------------------

float
computeLeastFloatSquaresEdgeLengths(const StrFloMatrix &orig_dm,  SequenceTree &tree){


  StrFloMatrix dm(orig_dm);
  const int numOriginalLeafs = dm.getSize();
  SequenceTree::NodeVector nodes;
  tree.recalcNodeIdsPostfixOrderAndAddInOrder(nodes);
  size_t nodeIdToRowIndex[nodes.size()];
  size_t rowIndexToNodeId[nodes.size()];
  str2int_hashmap name2Id((int)(nodes.size()*1.7));

  for(size_t i=0 ; i<nodes.size() ; i++)
    if(nodes[i]->isLeaf()){
      //PRINT(NAME(nodes[i]));PRINT(ID(nodes[i]));
      name2Id[NAME(nodes[i])] = ID(nodes[i]);
    }

  for(size_t row=0 ; row<dm.getSize() ; row++){
    str2int_hashmap::iterator f = name2Id.find(dm.getIdentifier(row));
    if(f==name2Id.end())
      USER_ERROR("name doesn't exist in tree: " << dm.getIdentifier(row));

    nodeIdToRowIndex[(*f).second] = row;
    rowIndexToNodeId[row] = (*f).second;
  }

  //the number of leafs below each node
  int numNodesBelow[nodes.size()];
  for(size_t i=0;i<nodes.size();i++)
    numNodesBelow[i]=1;

  //--------------------------------
  //BOTTOM UP TRAVERSAL IN TREE
  for(size_t i=0;i<nodes.size()-1;i++){
    if(nodes[i]->isLeaf())
      continue;

    //get the children and do the UNJ calculation to get the edge lengths
    SequenceTree::Node *parent = nodes[i];
    SequenceTree::Node *child1 = parent->getRightMostChild();
    SequenceTree::Node *child2 = child1->getLeftSibling();
    if(child2->getLeftSibling()!=NULL ){
      USER_ERROR("Have to be unrooted binary tree. Parent has " << parent->getNumChildren() << " children");
    }
    numNodesBelow[ID(parent)] = numNodesBelow[ID(child1)] + numNodesBelow[ID(child2)];
    //SEPARATOR();PRINT(NAME(child1));PRINT(NAME(child2));


    double sum = 0;
    for(size_t row=0;row<dm.getSize();row++){
      if(row==nodeIdToRowIndex[ID(child1)] ||
	 row==nodeIdToRowIndex[ID(child2)] )
	continue;
      sum += numNodesBelow[rowIndexToNodeId[row]]*(dm.getDistance(nodeIdToRowIndex[ID(child1)],row)-
						   dm.getDistance(nodeIdToRowIndex[ID(child2)],row));
    }

    if(!isfinite(sum)){
      USER_ERROR("Distance Matrix contains a non finite number: " << sum);
    }

    EDGE(child1) = 0.5*dm.getDistance(nodeIdToRowIndex[ID(child1)],
				      nodeIdToRowIndex[ID(child2)])
      + 1.0/(2*(numOriginalLeafs-numNodesBelow[ID(parent)]))*sum;
    EDGE(child2) = 0.5*dm.getDistance(nodeIdToRowIndex[ID(child1)],
				      nodeIdToRowIndex[ID(child2)])
      - 1.0/(2*(numOriginalLeafs-numNodesBelow[ID(parent)]))*sum;
    // PRINT(dm.getDistance(nodeIdToRowIndex[ID(child1)],nodeIdToRowIndex[ID(child2)]));
    //     PRINT((numOriginalLeafs-numNodesBelow[ID(parent)]));
    //     PRINT(sum);PRINT(EDGE(child1));PRINT( EDGE(child2));
    //PRINT(1/(2*(numOriginalLeafs-numNodesBelow[ID(parent)]))*sum);

    //swap child1 to last row
    int idOnLastRow = rowIndexToNodeId[dm.getSize()-1];
    if(idOnLastRow!=ID(child1)){
      int rowChild1 = nodeIdToRowIndex[ID(child1)];
      //PRINT(nodeIdToRowIndex[ID(child1)]);PRINT(dm.getSize());
      dm.swapRowToLast(nodeIdToRowIndex[ID(child1)]);
      nodeIdToRowIndex[idOnLastRow] = rowChild1;
      rowIndexToNodeId[rowChild1] = idOnLastRow;
      rowIndexToNodeId[dm.getSize()-1] = ID(child1);
      nodeIdToRowIndex[ID(child1)] = dm.getSize()-1;
    }
    //update distances to parent
    float w1 = (1.0*numNodesBelow[ID(child1)])/numNodesBelow[ID(parent)];
    float w2 = (1.0*numNodesBelow[ID(child2)])/numNodesBelow[ID(parent)];
    float distChild1Child2 = w1*EDGE(child1)+w2*EDGE(child2);

    //put parent on the row of child 2
    nodeIdToRowIndex[ID(parent)] = nodeIdToRowIndex[ID(child2)];
    rowIndexToNodeId[nodeIdToRowIndex[ID(parent)]] = ID(parent);
    int parentRow = nodeIdToRowIndex[ID(parent)];
    int child1Row = nodeIdToRowIndex[ID(child1)];
    int child2Row = nodeIdToRowIndex[ID(child2)];

    for(size_t row=0 ; row<dm.getSize()-1 ; row++){
      dm.setDistance(parentRow,row,
		     w1*dm.getDistance(child1Row,row)+
		     w2*dm.getDistance(child2Row,row)-
		     distChild1Child2);
    }

    dm.setDistance(nodeIdToRowIndex[ID(parent)],nodeIdToRowIndex[ID(parent)],0.0);
    //remove last row
    dm.removeLastRow();
  }

  //Take care of root
  SequenceTree::Node *root = nodes[nodes.size()-1];
  if(!root->isRoot() || root->getNumChildren()!=3){
    USER_ERROR("Have to be unrooted binary tree. Root has " << root->getNumChildren() << " children");
  }

  //  cout << dm << endl;
  SequenceTree::Node *c1 = root->getRightMostChild();
  SequenceTree::Node *c2 = c1->getLeftSibling();
  SequenceTree::Node *c3 = c2->getLeftSibling();

  int c1row = nodeIdToRowIndex[ID(c1)];
  int c2row = nodeIdToRowIndex[ID(c2)];
  int c3row = nodeIdToRowIndex[ID(c3)];

  EDGE(c1) = 0.5*(dm.getDistance(c1row,c2row) + dm.getDistance(c1row,c3row)-dm.getDistance(c2row,c3row));
  EDGE(c2) = 0.5*(dm.getDistance(c2row,c1row) + dm.getDistance(c2row,c3row)-dm.getDistance(c1row,c3row));
  EDGE(c3) = 0.5*(dm.getDistance(c3row,c2row) + dm.getDistance(c3row,c1row)-dm.getDistance(c2row,c1row));
  EDGE(root) = 0;

  //COMPUTE THE L2 SCORE
  StrFloMatrix treeM(tree.getNumLeafs());
  tree.tree2FloatdistanceMatrix(treeM);
  return computeFloatL2(treeM, orig_dm);
}


//--------------------------------------

float
computeFloatL2(const StrFloMatrix &A,  const StrFloMatrix &B){

  ASSERT_EQ(A.getSize(),B.getSize());

  str2int_hashmap name2Row((int)(A.getSize()*1.7));
  for(size_t i=0 ; i<A.getSize() ; i++){
    name2Row[A.getIdentifier(i)] = i;
  }

  float l2sum = 0;
  for(size_t i=0 ; i<B.getSize() ; i++){
    for(size_t j=i+1 ; j<B.getSize() ; j++ ){
      float b = B.getDistance(i,j);
      float a = A.getDistance(name2Row[B.getIdentifier(i)],
			       name2Row[B.getIdentifier(j)]);
      l2sum += (a-b)*(a-b);
    }
  }

  return l2sum;
}


// ----------- warka dang khatam sho------------
