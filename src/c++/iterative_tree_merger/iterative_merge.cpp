
#include "Tree.hpp"
#include "DistanceMatrix.hpp"

#include "FloatDistanceMatrix.hpp"
#include <string>

#include "arg_utils.h"
#include "log_utils.hpp"
#include <iostream>
#include <float.h>
#include <vector>
#include <algorithm>
#include "stl_utils.hpp"
#include "file_utils.hpp"
#include <fstream>
#include "NeighborJoining.hpp"
using namespace std;

//-------------------------------------------------------------
typedef DistanceMatrix< string, double, Data_init<string>, Data_printOn<string>, Data_init<double>, Data_printOn<double> > DM;
typedef Tree<string, Data_init<string>, Data_printOn<string> > StringTree;
typedef Tree<string_int, Data_init<string_int>, Data_printOn<string_int> > StringIntTree;
typedef DistanceMatrix< StringIntTree::Node *, double, Data_init<StringIntTree::Node *>, Data_printOn<StringIntTree::Node *>, Data_init<double>, Data_printOn<double> > StringIntTree_DM;
typedef DistanceMatrix< StringTree::Node *, double, Data_init<StringTree::Node *>, Data_printOn<StringTree::Node *>, Data_init<double>, Data_printOn<double> > StringTree_DM;

struct node_dist{
  StringIntTree::Node *n1;
  StringIntTree::Node *n2;
  double dist;
};

void
NJ_MERGE_TREES(StringIntTree *t1, StringIntTree *t2, StringIntTree_DM &dm);

bool
operator < (const node_dist &n1, const node_dist &n2){
  return n1.dist < n2.dist;
}


typedef __gnu_cxx::hash_set<StringIntTree::Node *, objhash_ptr, objeq_ptr> node_set;
typedef __gnu_cxx::hash_map<StringIntTree::Node *, int, objhash_ptr, objeq_ptr> node_int_map;


//---------------------------------
StringTree
computeIterativeMergeTree( DM &orig_dm ){

  int n = orig_dm.getSize();
  //PRINT(orig_dm);
  //create a tree for each string
  vector<StringIntTree *> trees;
  StringIntTree::NodeVector nodes;
  for ( int i=0 ; i<n ; i++ ){
    string_int strint;
    strint.s = orig_dm.getIdentifier(i);
    strint.i = i;
    StringIntTree *t = new StringIntTree(strint);
    trees.push_back(t);
    nodes.push_back(t->getRoot());
  }

  //copy the distances 
  vector<node_dist> node_dists;
  StringIntTree_DM dm(n);
  for ( int i=0 ; i<n ; i++ ){
    dm.setIdentifier(i,nodes[i]);
    dm.setDistance(i,i,0);
  }
  for ( int i=0 ; i<n ; i++ ){
    for ( int j=i+1 ; j<n ; j++ ){
      dm.setDistance(i,j,orig_dm.getDistance(i,j));
      node_dist nd;
      nd.n1 = nodes[i];
      nd.n2 = nodes[j];
      nd.dist = orig_dm.getDistance(i,j);
      node_dists.push_back(nd);
    } 
  }

  //----
  //sort distances
  sort(node_dists.begin(),node_dists.end());
  size_t dist_index =0;
  //-----
  //THE MERGER LOOP
  for ( int i=0 ; i<n-1 ; i++ ){
    //PRINT(i);
    //take the minimum leaf to leaf distance with leafs of different trees
    node_dist nd;
    do{
      nd = node_dists[dist_index];
      // PRINT(nd.n1);PRINT(nd.n2);
      dist_index++;
    }while(nd.n1->getTree() == nd.n2->getTree());
    // PRINT(i);
    //     PRINT(nd.n1->getTree()->getNumNodes());
    //     PRINT(nd.n2->getTree()->getNumNodes());
    //     PRINT(nd.n1->getTree());
    //     PRINT(nd.n2->getTree());
    
    NJ_MERGE_TREES(nd.n1->getTree(), nd.n2->getTree(), dm);
    //nd.n1->getTree()->drawTree(cout);
    //PRINT(nd.n1->getTree());
    //PRINT(nd.n2->getTree());
    
  }

  nodes[0]->getTree()->reRootAtHighDegreeNode();
  ASSERT_EQ((int) orig_dm.getSize(), (int) nodes[0]->getTree()->getNumLeafs() );
  //CREATE THE STRING TREE
  StringIntTree *tree = nodes[0]->getTree();
  tree->recalcNodeStructure();
  StringIntTree::NodeVector leafs;
  tree->addLeafs(leafs);
    
  string dummystr;
  StringTree resultTree(*tree,dummystr);
  StringTree::NodeVector resultleafs;
  resultTree.addLeafs(resultleafs);

  for ( size_t i=0 ; i<leafs.size() ; i++ ){
    resultleafs[i]->data = leafs[i]->data.s;
  }

  //FREE THE MEMORY
  for ( size_t i=0 ; i<trees.size() ;i++ )
    delete trees[i];

  //   //DEBUGGING
  //   //----- DO NJ
  //   string s ="";
  //   StringTree_DM njdm(orig_dm.getSize());
  //   StringTree resultTree = StringTree(s);
  //   for ( size_t i=0 ; i<orig_dm.getSize() ; i++ ){
  //     StringTree::Node *n = resultTree.getRoot()->addChild(orig_dm.getIdentifier(i));
  //     njdm.setIdentifier(i,n);
  //     for ( size_t j = i ; j<orig_dm.getSize() ; j++ )
  //       njdm.setDistance(i,j,orig_dm.getDistance(i,j));
  //   }
  //   computeNeighborJoiningTree(njdm,s);
  //-----
  return resultTree;
}


//-------------------------
void
NJ_MERGE_TREES(StringIntTree *t1, StringIntTree *t2, StringIntTree_DM &big_dm){
  
  vector<StringIntTree::Node *> leafs;
  if ( t1->getNumNodes() == 1 )
    leafs.push_back(t1->getRoot());
  else
    t1->addNodesWithDegree(leafs,1);

  if ( t2->getNumNodes() == 1 )
    leafs.push_back(t2->getRoot());
  else
    t2->addNodesWithDegree(leafs,1);


  if ( leafs.size() == 2 ){
    t1->joinTreeAtRoot(*t2);
    return;
  }
  //  PRINT(leafs.size());

  //-----
  //copy distances between leafs
  StringIntTree_DM dm(leafs.size());
  for ( size_t i=0 ; i<leafs.size() ; i++ ){
    dm.setIdentifier(i,leafs[i]);
    dm.setDistance(i,i,0);
  }

  for ( size_t i=0 ; i<leafs.size() ; i++ ){
    for ( size_t j=0 ; j<leafs.size() ; j++ ){
      dm.setDistance(i,j,big_dm.getDistance(leafs[i]->data.i,leafs[j]->data.i));
    }
  }
  //  PRINT(dm);
  //------
  //START NJ

  //compute the row sums
  double rowSums[dm.getSize()];
  for ( size_t row = 0 ; row < dm.getSize() ; row++ ){
    double sum = 0;
    size_t i = 0;
    for (  ; i < dm.getSize() ; i++ )
      sum += dm.getDistance(row,i);
    rowSums[row] = sum;
  }


  //----------------
  // NJ ITERATION
  string_int si;//default string_int
  si.s="";si.i=-1;
  node_set usednodes;
  while ( dm.getSize() > 3 ) {
    //find the minimal value
    double minVal = FLT_MAX;
    size_t mini = 1000000;
    size_t minj = 1000000;
    for ( size_t i = 0 ; i < dm.getSize() ; i++ ){
      for ( size_t j = i+1 ; j < dm.getSize() ; j++ ){
        double newVal = (dm.getSize() - 2.0)*dm.getDistance(i,j) - rowSums[i] - rowSums[j];
        if ( newVal < minVal ){

          //if the nodes are in the same tree check that they are neighbors
          StringIntTree::Node *n1 = dm.getIdentifier(i);
          StringIntTree::Node *n2 = dm.getIdentifier(j);
          if ( n1->getTree() == n2->getTree() ){
            StringIntTree::Node *n1p = n1->getParent();
            StringIntTree::Node *n2p = n2->getParent();
            if ( n1p!=NULL && 
                 (n1p==n2 || n1p==n2p || n1p->getParent() == n2 ) )
              goto ok_neighbors;
            if ( n2p!=NULL && 
                 (n2p==n1 || n2p->getParent() == n1 ) )
              goto ok_neighbors;
            //not ok breaks neighbor condition
            continue;
          }
        ok_neighbors:
          minVal = newVal;
          mini = i;
          minj = j;
        }
      }
    }

    //make sure that minj is the last row in the matrix
    if ( minj != dm.getSize() -1 ){
      dm.swapRowToLast(mini);
      double tmp = rowSums[dm.getSize()-1];
      rowSums[dm.getSize()-1] = rowSums[mini];
      rowSums[mini] = tmp;
      mini = minj;
      minj = dm.getSize()-1;
    }

    //join the tree nodes
    StringIntTree::Node *n1 = dm.getIdentifier(mini);
    StringIntTree::Node *n2 = dm.getIdentifier(minj);
    //CASE 1 if different trees
    if ( n1->getTree() != n2->getTree() ){
      //If degree 0 then just connect to new node
      //if degree 1 then infer new node on path and connect
      //if degree 2 (then root node) just connect.
      //if degree 3 then one of the neighbors have not been used.
      //connect the trees on the edge to the unused node 
      bool crash = true;
      
      StringIntTree::Node *connect1 = NULL;
      if ( n1->getDegree() == 0 || n1->getDegree()==2){
        crash = false;
        connect1 = n1;
      }
      else if ( n1->getDegree() == 1 ){
        if ( n1->getParent() == NULL )
          n1->getTree()->reRootAt(n1->getRightMostChild());
	
        connect1 = n1->getTree()->insertNodeOnPathToParent(n1, si);
      }
      else if ( n1->getDegree() == 3 ){
        n1->getTree()->reRootAt(n1);
        connect1 = n1->getLeftMostChild();
        while ( connect1 != NULL ){
          if ( usednodes.find(connect1) == usednodes.end() )
            break;
          connect1 = connect1->getRightSibling();
        }
        n1->getTree()->reRootAt(connect1);
        connect1 = n1->getTree()->insertNodeOnPathToParent(n1, si);
      }
      assert(connect1 != NULL );
      //--
      StringIntTree::Node *connect2 = NULL;
      if ( n2->getDegree() == 0 || n2->getDegree()==2){
        crash = false;
        connect2 = n2;
      }
      else if ( n2->getDegree() == 1 ){
        if ( n2->getParent() == NULL )
          n2->getTree()->reRootAt(n2->getRightMostChild());
        connect2 = n2->getTree()->insertNodeOnPathToParent(n2, si);
      }
      else if ( n2->getDegree() == 3 ){
        n2->getTree()->reRootAt(n2);
        connect2 = n2->getLeftMostChild();
        while ( connect2 != NULL ){
          if ( usednodes.find(connect2) == usednodes.end() )
            break;
          connect2 = connect2->getRightSibling();
        }
        n2->getTree()->reRootAt(connect2);

        connect2 = n2->getTree()->insertNodeOnPathToParent(n2, si);
      }
      assert(connect2 != NULL );
      //if( crash ) cout << "CRASH" << endl;
      //--------
      //FOUND MERGE POINT
      //collapse the connect nodes so there are 16 children of the root
      connect1->getTree()->reRootAt(connect1);
      connect2->getTree()->reRootAt(connect2);
      connect1->getTree()->collapseChildren(connect1);
      connect1->getTree()->collapseChildren(connect1);
      connect2->getTree()->collapseChildren(connect2);
      connect2->getTree()->collapseChildren(connect2);
      StringIntTree *tree;
      if(connect1->getDegree()>=2){
        connect1->getTree()->joinTreeAtRoot(*(connect2->getTree()));
        connect2->getTree()->collapse(connect2);
        tree = connect1->getTree();
      }
      else {
        connect2->getTree()->joinTreeAtRoot(*(connect1->getTree()));
        connect1->getTree()->collapse(connect1);
        tree = connect2->getTree();
      }
      //now connect1 or connect2 is the root and has <=16 children
      //---------
      //DO NJ on the 16-star after collapse
      //compute the distances for the children of the root
      //copy distances between leafs
      dm.resize(leafs.size());
      node_int_map rowMap((int)(leafs.size()*1.3));
      for ( size_t i=0 ; i<leafs.size() ; i++ ){
        dm.setIdentifier(i,leafs[i]);
        rowMap[leafs[i]] = i;
        dm.setDistance(i,i,0);
      }
      
      for ( size_t i=0 ; i<leafs.size() ; i++ ){
        for ( size_t j=0 ; j<leafs.size() ; j++ ){
          dm.setDistance(i,j,big_dm.getDistance(leafs[i]->data.i,leafs[j]->data.i));
        }
      }

      //do nj updates for the nodes in the tree.
      StringIntTree::NodeVector postfix;
      tree->addNodesInPostfixOrder(postfix);
      for ( size_t i = 0 ; i < postfix.size()-1 ; i++ ){
        if ( postfix[i]->isLeaf() ) continue;
        StringIntTree::Node *rightChild = postfix[i]->getRightMostChild();
        StringIntTree::Node *leftChild = rightChild->getLeftSibling();
        ASSERT_EQ(leftChild,postfix[i]->getLeftMostChild());//should only have to children
        size_t rowR = rowMap[rightChild];
        size_t rowL = rowMap[leftChild];
        if ( rowL == dm.getSize()-1 )
          rowL = rowR;//warning... not changing the identifier
        else{
          rowMap[dm.getIdentifier(dm.getSize()-1)] = rowR;
          dm.swapRowToLast(rowR);
        }
        rowR = dm.getSize()-1;
        for ( size_t j = 0 ; j < rowR ; j++ )
          dm.setDistance(rowL,j,(dm.getDistance(rowL,j)+dm.getDistance(rowR,j))/2.0);
        dm.setDistance(rowL,rowL,0);
        dm.setIdentifier(rowL,postfix[i]);
        dm.removeLastRow();
        rowMap[postfix[i]] = rowL;
      }
      //compute NJ for the 16 star
      
      computeNeighborJoiningTree(dm,si);
      return;//FOUND MERGE POINT
    }
    //---------------------------------------
    //CASE 2 the same tree
    //either replace by old node or finnished the tree and has to infer new node on the path.
    StringIntTree::Node *n1p = n1->getParent();
    StringIntTree::Node *n2p = n2->getParent();
    StringIntTree::Node *newparent;
    if ( n1p!=NULL && n1p == n2  ){
      //FINNISHED tree
      newparent = n1->getTree()->insertNodeOnPathToParent(n1, si);
    }
    else if ( n2p!=NULL && n2p == n1 ){
      //FINNISHED tree
      newparent = n2->getTree()->insertNodeOnPathToParent(n2, si);
    }
    else if ( n1p == n2p )
      newparent = n1p;
    else if ( n1p!=NULL && n1p->getParent() == n2 )
      newparent = n1p;
    else 
      newparent = n2p;

    usednodes.insert(n1);
    usednodes.insert(n2);

    dm.setIdentifier(mini, newparent);
    
    // UPDATE DISTANCES
    for ( size_t i = 0 ; i < dm.getSize()-1 ; i++ ){//skip last row
      double dist2iandj = dm.getDistance(mini,i) + dm.getDistance(minj,i);
      // regular nj update function:
      dm.setDistance(mini,i, dist2iandj * 0.5); 
    
      //update rowsums
      rowSums[i] = rowSums[i] - dist2iandj + dm.getDistance(mini,i);
    }

    //remove the last row of the matrix
    dm.removeLastRow();
    
    //recompute the row sum for the parent
    dm.setDistance(mini,mini,0);
    double sum = 0;
    for ( size_t i = 0 ; i < dm.getSize() ; i++ )
      sum += dm.getDistance(mini,i);
    rowSums[mini] = sum;
  }
  // END ITERATION
  //--------------

  //merge the three subtrees
  //PRINT(dm.getSize());
  StringIntTree::Node *n1 = dm.getIdentifier(0);
  StringIntTree::Node *n2 = dm.getIdentifier(1);
  StringIntTree::Node *n3 = dm.getIdentifier(2);
  StringIntTree::Node *newparent;
  if ( n1->getTree() == n2->getTree() ){
    if ( n1->getParent() == n2 ){
      newparent = n1->getTree()->insertNodeOnPathToParent(n1, si);
    }
    else if ( n2->getParent() == n1 ){
      newparent = n2->getTree()->insertNodeOnPathToParent(n2, si);
    }
    else{
      PROG_ERROR("shouldn't come here");
    }
    n1->getTree()->joinTreeAtNodes(newparent,n3); 
  }
  else if ( n1->getTree() == n3->getTree() ){
    if ( n1->getParent() == n3 ){
      newparent = n1->getTree()->insertNodeOnPathToParent(n1, si);
    }
    else if ( n3->getParent() == n1 ){
      newparent = n3->getTree()->insertNodeOnPathToParent(n3, si);
    }
    else{
      PROG_ERROR("shouldn't come here");
    }
    n1->getTree()->joinTreeAtNodes(newparent,n2); 
  }
  else {//
    ASSERT_EQ( n2->getTree(), n3->getTree() );
    if ( n2->getParent() == n3 ){
      newparent = n2->getTree()->insertNodeOnPathToParent(n2, si);
    }
    else if ( n3->getParent() == n2 ){
      newparent = n3->getTree()->insertNodeOnPathToParent(n3, si);
    }
    else{
      PROG_ERROR("shouldn't come here");
    }
    //PRINT(n2->getTree());
    n2->getTree()->joinTreeAtNodes(newparent,n1); 
  }
}







void
print_usage(){
  cout << "usage: interative_merge -i dmatrix -m 20"  <<endl;
}

int
main(int argc, char **argv){

  char *dfile = GET_OPTION_VAL("-i");
  if ( dfile == NULL ) print_usage();
  ifstream in;
  open_read_stream(dfile,in);

  int m = 1;
  char *mstr = GET_OPTION_VAL("-m");
  if ( mstr != NULL)
    m = atoi(mstr);

  for ( int i=0 ; i<m ; i++ ){
    DM dm(in);

    StringTree tree = computeIterativeMergeTree(dm);

    cout << tree << endl;
  }
  in.close();
}



