
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
#include <cmath>
#include <vector>
#include "fastphylo/core/SequenceTree.hpp"
#include "fastphylo/dna/LeastSquaresFit.hpp"

using namespace std;

#define ID(node) (node)->getNodeId()

// computeLeastSquaresEdgeLengths()/computeLeastFloatSquaresEdgeLengths()
// below run the identical UNJ edge-length-fitting algorithm, one on
// StrDblMatrix/double and one on StrFloMatrix/float - the double/float
// distinction only matters for the matrix type itself, the w1/w2/
// distChild1Child2 locals (Real below), and which of tree2distanceMatrix/
// tree2FloatdistanceMatrix and computeL2/computeFloatL2 to call (different
// names, not overloads, so passed in as function pointers). Shared here
// as a template so the actual algorithm exists in one place instead of
// two 135-line copies.
namespace {

// Maps each leaf's node id to dm's row index and back, so the bottom-up
// traversal below can go by node id instead of re-searching dm each time.
template <typename Matrix>
void buildRowIndexMaps(const Matrix &dm, const SequenceTree::NodeVector &nodes,
                        std::vector<size_t> &nodeIdToRowIndex, std::vector<size_t> &rowIndexToNodeId){
  str2int_hashmap name2Id(static_cast<int>(nodes.size()*1.7));
  for(size_t i=0 ; i<nodes.size() ; i++) {
    if(nodes[i]->isLeaf()){
      //PRINT(NAME(nodes[i]));PRINT(ID(nodes[i]));
      name2Id[NAME(nodes[i])] = ID(nodes[i]);
    }
}

  for(size_t row=0 ; row<dm.getSize() ; row++){
    auto f = name2Id.find(dm.getIdentifier(row));
    if(f==name2Id.end())
      USER_ERROR("name doesn't exist in tree: " << dm.getIdentifier(row));

    nodeIdToRowIndex[(*f).second] = row;
    rowIndexToNodeId[row] = (*f).second;
  }
}

// The UNJ weighted-sum term used by both of EDGE(child1)/EDGE(child2)'s
// formulas below - every other row's distance to child1 minus its
// distance to child2, weighted by how many leaves are below that row.
template <typename Matrix>
double computeWeightedNeighborSum(const Matrix &dm, const std::vector<int> &numNodesBelow,
                                   const std::vector<size_t> &rowIndexToNodeId,
                                   size_t child1Row, size_t child2Row){
  double sum = 0;
  for(size_t row=0;row<dm.getSize();row++){
    if(row==child1Row || row==child2Row) {
      continue;
}
    sum += numNodesBelow[rowIndexToNodeId[row]]*(dm.getDistance(child1Row,row)-
						 dm.getDistance(child2Row,row));
  }
  return sum;
}

// Moves child1's row to the last row of dm (so removeLastRow() below
// drops the right one), updating the row<->id maps to match.
template <typename Matrix>
void swapChildToLastRow(Matrix &dm, std::vector<size_t> &nodeIdToRowIndex,
                         std::vector<size_t> &rowIndexToNodeId, int child1Id){
  int idOnLastRow = rowIndexToNodeId[dm.getSize()-1];
  if(idOnLastRow!=child1Id){
    int rowChild1 = nodeIdToRowIndex[child1Id];
    //PRINT(nodeIdToRowIndex[child1Id]);PRINT(dm.getSize());
    dm.swapRowToLast(nodeIdToRowIndex[child1Id]);
    nodeIdToRowIndex[idOnLastRow] = rowChild1;
    rowIndexToNodeId[rowChild1] = idOnLastRow;
    rowIndexToNodeId[dm.getSize()-1] = child1Id;
    nodeIdToRowIndex[child1Id] = dm.getSize()-1;
  }
}

// Replaces child1's/child2's rows in dm with parent's UNJ-weighted
// combined distance to every other (not-yet-collapsed) row.
template <typename Matrix, typename Real>
void updateDistancesToParent(Matrix &dm, size_t parentRow, size_t child1Row, size_t child2Row,
                              Real w1, Real w2, Real distChild1Child2){
  for(size_t row=0 ; row<dm.getSize()-1 ; row++){
    dm.setDistance(parentRow,row,
		   (w1*dm.getDistance(child1Row,row))+
		   (w2*dm.getDistance(child2Row,row))-
		   distChild1Child2);
  }
}

template <typename Matrix, typename Real>
Real computeLeastSquaresEdgeLengthsImpl(const Matrix &orig_dm, SequenceTree &tree,
                                         void (SequenceTree::*treeToMatrix)(Matrix &),
                                         Real (*computeL2Fn)(const Matrix &, const Matrix &)){

  Matrix dm(orig_dm);
  const int numOriginalLeafs = dm.getSize();
  SequenceTree::NodeVector nodes;
  tree.recalcNodeIdsPostfixOrderAndAddInOrder(nodes);
  std::vector<size_t> nodeIdToRowIndex(nodes.size());
  std::vector<size_t> rowIndexToNodeId(nodes.size());
  buildRowIndexMaps(dm, nodes, nodeIdToRowIndex, rowIndexToNodeId);

  //the number of leafs below each node
  std::vector<int> numNodesBelow(nodes.size());
  for(size_t i=0;i<nodes.size();i++) {
    numNodesBelow[i]=1;
}

  //--------------------------------
  //BOTTOM UP TRAVERSAL IN TREE
  for(size_t i=0;i<nodes.size()-1;i++){
    if(nodes[i]->isLeaf()) {
      continue;
}

    //get the children and do the UNJ calculation to get the edge lengths
    SequenceTree::Node *parent = nodes[i];
    SequenceTree::Node *child1 = parent->getRightMostChild();
    SequenceTree::Node *child2 = child1->getLeftSibling();
    if(child2->getLeftSibling()!=nullptr ){
      USER_ERROR("Have to be unrooted binary tree. Parent has " << parent->getNumChildren() << " children");
    }
    numNodesBelow[ID(parent)] = numNodesBelow[ID(child1)] + numNodesBelow[ID(child2)];
    //SEPARATOR();PRINT(NAME(child1));PRINT(NAME(child2));


    double sum = computeWeightedNeighborSum(dm, numNodesBelow, rowIndexToNodeId,
                                             nodeIdToRowIndex[ID(child1)], nodeIdToRowIndex[ID(child2)]);

    if(!isfinite(sum)){
      USER_ERROR("Distance Matrix contains a non finite number: " << sum);
    }

    EDGE(child1) = (0.5*dm.getDistance(nodeIdToRowIndex[ID(child1)],
				      nodeIdToRowIndex[ID(child2)]))
      + (1.0/(2*(numOriginalLeafs-numNodesBelow[ID(parent)]))*sum);
    EDGE(child2) = (0.5*dm.getDistance(nodeIdToRowIndex[ID(child1)],
				      nodeIdToRowIndex[ID(child2)]))
      - (1.0/(2*(numOriginalLeafs-numNodesBelow[ID(parent)]))*sum);
    // PRINT(dm.getDistance(nodeIdToRowIndex[ID(child1)],nodeIdToRowIndex[ID(child2)]));
    //     PRINT((numOriginalLeafs-numNodesBelow[ID(parent)]));
    //     PRINT(sum);PRINT(EDGE(child1));PRINT( EDGE(child2));
    //PRINT(1/(2*(numOriginalLeafs-numNodesBelow[ID(parent)]))*sum);

    //swap child1 to last row
    swapChildToLastRow(dm, nodeIdToRowIndex, rowIndexToNodeId, ID(child1));
    //update distances to parent
    Real w1 = (1.0*numNodesBelow[ID(child1)])/numNodesBelow[ID(parent)];
    Real w2 = (1.0*numNodesBelow[ID(child2)])/numNodesBelow[ID(parent)];
    Real distChild1Child2 = (w1*EDGE(child1))+(w2*EDGE(child2));

    //put parent on the row of child 2
    nodeIdToRowIndex[ID(parent)] = nodeIdToRowIndex[ID(child2)];
    rowIndexToNodeId[nodeIdToRowIndex[ID(parent)]] = ID(parent);
    int parentRow = nodeIdToRowIndex[ID(parent)];
    int child1Row = nodeIdToRowIndex[ID(child1)];
    int child2Row = nodeIdToRowIndex[ID(child2)];

    updateDistancesToParent(dm, parentRow, child1Row, child2Row, w1, w2, distChild1Child2);

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
  Matrix treeM(tree.getNumLeafs());
  (tree.*treeToMatrix)(treeM);
  return computeL2Fn(treeM, orig_dm);
}
} // namespace

double
computeLeastSquaresEdgeLengths(const StrDblMatrix &orig_dm,  SequenceTree &tree){
  return computeLeastSquaresEdgeLengthsImpl<StrDblMatrix, double>(orig_dm, tree, &SequenceTree::tree2distanceMatrix, &computeL2);
}


//--------------------------------------

double
computeL2(const StrDblMatrix &A,  const StrDblMatrix &B){
  
  ASSERT_EQ(A.getSize(),B.getSize());
  
  str2int_hashmap name2Row(static_cast<int>(A.getSize()*1.7));
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
  return computeLeastSquaresEdgeLengthsImpl<StrFloMatrix, float>(orig_dm, tree, &SequenceTree::tree2FloatdistanceMatrix, &computeFloatL2);
}


//--------------------------------------

float
computeFloatL2(const StrFloMatrix &A,  const StrFloMatrix &B){

  ASSERT_EQ(A.getSize(),B.getSize());

  str2int_hashmap name2Row(static_cast<int>(A.getSize()*1.7));
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
