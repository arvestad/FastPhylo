//--------------------------------------------------
//                                        
// File: FNJ_collapse_method.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: FNJ_collapse_method.cpp,v 1.5 2006/08/08 13:50:51 isaac Exp $                                 
//
//-------------------------------------------------

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
#include <math.h>

using namespace std;


//---------------------------
// The Tree
// Each node in the tree has a data containing a string and an int.
// The string represents the name of the taxa
// and the int represents the rowId in the distance matrix
typedef Tree<string_int, Data_init<string_int>, Data_printOn<string_int> > StringIntTree;

typedef Tree<int_double, Data_init<int_double>, Data_printOn<int_double> >  IntDoubleTree;

//------------------
// The distance matrix
//
// Each row has a tree node as an identifier.  The matrix is of
// dimension (2n)*(2n).  OBSERVE: All nodes in the tree are part of the
// matrix; in regular NJ only the children of the root are part of the
// tree.
typedef DistanceMatrix< StringIntTree::Node *, double, Data_init<StringIntTree::Node *>, Data_printOn<StringIntTree::Node *>, Data_init<double>, Data_printOn<double> > StringIntTree_DM;

//The input distance matrix
typedef DistanceMatrix< string, double, Data_init<string>, Data_printOn<string>, Data_init<double>, Data_printOn<double> > String_DM;

typedef DistanceMatrix< int, double, Data_init<int>, Data_printOn<int>, Data_init<double>, Data_printOn<double> > Int_DM;

//--------------------------------------
StringIntTree tree;
StringIntTree::Node *root;
StringIntTree_DM dm(2);
StringIntTree_DM varianceMatrix(2);
vector<double> activeRowSums;
vector<double> varianceRowSums;
vector<int> vacantRows;
int numChildrenOfRoot;
int numLeafs;
int sequenceLength;//used by weighbor
//-----------------------------------------------------------

int COLLAPSE_LEVELS=3;
char SUBTREE_METHOD='B';


void
run_subtree_WeighborExternally_on_children(StringIntTree::Node *parent){

  //build distance matrix for children
  int numC = parent->getNumChildren();
  Int_DM tmp_dm(numC);
  StringIntTree::Node *child1 = parent->getRightMostChild();
  for( int child_i=0 ; child_i<numC ; child_i++,child1=child1->getLeftSibling() ){
    tmp_dm.setIdentifier(child_i,child1->data.i);
    StringIntTree::Node *child2 = child1->getLeftSibling();
    for( int child_j=child_i+1 ; child_j<numC ; child_j++,child2=child2->getLeftSibling() )
      tmp_dm.setDistance(child_i,child_j,dm.getDistance(child1->data.i,child2->data.i));
  }

  ofstream dmfile;
  char *dmfilename = "xxMATRIXxx.dist";
  open_write_stream(dmfilename,dmfile);
  dmfile << tmp_dm << endl;
  dmfile.close();

  //call weighbor
  char *wresult = "weighTree";
  string wstr = string("~/10giga_volume/Weighbor/weighbor -L ") +sequenceLength + " -i " + string(dmfilename) +"  > " + wresult;
  system(wstr.c_str());

  //read the tree
  ifstream treefile;
  open_read_stream(wresult,treefile);
  IntDoubleTree wtree(treefile);
  treefile.close();
  //  wtree.drawTree(cout);

  //reconstruct the children of parent as described by weighbor
  IntDoubleTree::NodeVector nodes;
  wtree.addNodesInPostfixOrder(nodes);
  
  for( size_t i=0 ; i<nodes.size() ; i++ ){
    if ( nodes[i]->getParent()==NULL || nodes[i]->getParent()->getParent()==NULL)
      continue;
    //get sibling
    if ( nodes[i]->getLeftSibling()==NULL && nodes[i]->getRightSibling()!=NULL ){
      StringIntTree::Node *pchild1 = dm.getIdentifier(nodes[i]->data.i);
      StringIntTree::Node *pchild2 = dm.getIdentifier(nodes[i]->getRightSibling()->data.i);
      //PRINT(nodes[i]->getNodeId());PRINT(nodes[i]->getRightSibling()->getNodeId());PRINT(pchild1->data.i);PRINT(pchild2->data.i);
      string_int si;
      si.i = vacantRows.back();
      vacantRows.pop_back();
      StringIntTree::Node *inner_parent = tree.detachFromParentAndAddAsSiblings(pchild1,pchild2,si);
      dm.setIdentifier(si.i,inner_parent);
      nodes[i]->getParent()->data.i = si.i;
      
      //      tree.drawSubtree(cout,parent);
    }
  }
  //  tree.drawSubtree(cout,parent);
}

void
collapse_children_levels(StringIntTree::Node *parent, int numLevels){
  //collapse the children of parent to 4 levels and remove the
  //collapsed nodes from the distance matrix
  StringIntTree::Node *nextChild;
  for( int level=0 ; level<numLevels ; level++ ){
    StringIntTree::Node *child = parent->getRightMostChild();
    while( child!=NULL ){
      nextChild = child->getLeftSibling();
      if ( !child->isLeaf() ){
        dm.setIdentifier(child->data.i,NULL);
        vacantRows.push_back(child->data.i);
        tree.collapse(child);
      }
      child = nextChild;
    }
  }

}

void
run_regular_NJ_on_children(StringIntTree::Node *parent){
  //---
  //run NJ on the children of parent until parent only has two children
  int min_i = -1;
  int min_j = -1;
  int numChildren = parent->getNumChildren();
 
  for( ; numChildren > 2; numChildren--){
    int numIter=0;
    double minVal =  FLT_MAX;
    StringIntTree::Node *child1 = parent->getRightMostChild();
    StringIntTree::Node *child2;
    for( ; child1!=NULL ; child1=child1->getLeftSibling() ){
      child2 = child1->getLeftSibling();
      //child2 = parent->getRightMostChild();
      for( ; child2!=NULL ; child2=child2->getLeftSibling() ){
        if ( child2==child1 ) continue;
        numIter++;
        double newVal = (numChildren+(numChildrenOfRoot-1) - 2.0)*dm.getDistance(child1->data.i,child2->data.i) - activeRowSums[child1->data.i] - activeRowSums[child2->data.i];
        if( newVal<minVal ){
          minVal = newVal;
          min_i = child1->data.i;
          min_j = child2->data.i;
        }
      }
    }
    //join the min pair into a parent node
    child1 = dm.getIdentifier(min_i);
    child2 = dm.getIdentifier(min_j);
    string_int si;
    si.i = vacantRows.back();
    vacantRows.pop_back();
    //LINE();
    StringIntTree::Node *inner_parent = tree.detachFromParentAndAddAsSiblings(child1, child2, si);
    dm.setIdentifier(si.i,inner_parent);

    //UPDATE ALL DISTANCES to inner_parent
    //compute lambda BIONJ way
    double r_2=(numChildren+(numChildrenOfRoot-1) - 2.0);
    double lambda = 0.5 + (varianceRowSums[min_j] - varianceRowSums[min_i])/(2.0*r_2*varianceMatrix.getDistance(min_i,min_j));
    if ( lambda > 1 ) lambda = 1.0;
    else if ( lambda < 0 ) lambda = 0;
    else if ( varianceMatrix.getDistance(min_i,min_j) == 0 )//i.e. lambda==NaN
      lambda = 0.5;
    
    double d1u = 0.5*(dm.getDistance(min_i,min_j) + (activeRowSums[min_i]-activeRowSums[min_j])/r_2);
    double d2u = 0.5*(dm.getDistance(min_i,min_j) + (activeRowSums[min_j]-activeRowSums[min_i])/r_2);

    for( int i=0 ; i<2*numLeafs ; i++ ){
      if ( dm.getIdentifier(i) == NULL ) continue;
      double dui = lambda*dm.getDistance(min_i,i) + (1.0-lambda)*dm.getDistance(min_j,i) - lambda*d1u - (1.0-lambda)*d2u;
      dm.setDistance(inner_parent->data.i,i, dui);
      double vui = lambda*varianceMatrix.getDistance(min_i,i) + (1-lambda)*varianceMatrix.getDistance(min_j,i) - lambda*(1-lambda)*varianceMatrix.getDistance(min_i,min_j);
      varianceMatrix.setDistance(inner_parent->data.i,i, vui);
    }
      
    //update the rowsums for the children of parent
    StringIntTree::Node *child = parent->getRightMostChild();
    for( ; child!=NULL ; child=child->getLeftSibling() ){
      if ( child==inner_parent ) continue;
      double dist2iandj = dm.getDistance(min_i,child->data.i) + dm.getDistance(min_j,child->data.i);
      //update rowsums
      activeRowSums[child->data.i] += (dm.getDistance(inner_parent->data.i,child->data.i)-dist2iandj);
      
      double variance2iandj = varianceMatrix.getDistance(min_i,child->data.i) + varianceMatrix.getDistance(min_j,child->data.i);
      varianceRowSums[child->data.i] += (varianceMatrix.getDistance(inner_parent->data.i,child->data.i) - variance2iandj);
    }
      
    dm.setDistance(inner_parent->data.i,inner_parent->data.i,0);
    varianceMatrix.setDistance(inner_parent->data.i,inner_parent->data.i,0);

    //---
    //compute the rowsum for the new node
    double rowSum = 0;
    double varianceSum = 0;
    //compute the sum to the distance to the children of parent
    child = parent->getRightMostChild();
    for( ; child!=NULL ; child=child->getLeftSibling() ){
      if ( child==inner_parent ) continue;
      rowSum += dm.getDistance(child->data.i, inner_parent->data.i);
      varianceSum += varianceMatrix.getDistance(child->data.i, inner_parent->data.i);
    }
    //compute the sum to the children of root
    child = root->getRightMostChild();
    for( ; child!=NULL ; child=child->getLeftSibling() ){
      if ( child==parent ) continue;
      rowSum += dm.getDistance(child->data.i, inner_parent->data.i);
      varianceSum += varianceMatrix.getDistance(child->data.i, inner_parent->data.i);
    }
    activeRowSums[inner_parent->data.i] = rowSum;
    varianceRowSums[inner_parent->data.i] = varianceSum;
    
    //----

  }//end inner NJ loop
  //----
}




void
run_subtree_NJ_on_children(StringIntTree::Node *parent){
  //---
  //run NJ on the children of parent until parent only has two children
  int min_i =-1;
  int min_j =-1;
  int numChildren = parent->getNumChildren();
  for( ; numChildren > 3; numChildren--){
    double minVal =  FLT_MAX;
    StringIntTree::Node *child1 = parent->getRightMostChild();
    StringIntTree::Node *child2;
    for( ; child1!=NULL ; child1=child1->getLeftSibling() ){
      child2 = child1->getLeftSibling();
      //child2 = parent->getRightMostChild();
      for( ; child2!=NULL ; child2=child2->getLeftSibling() ){
        if ( child2==child1 ) continue;
        double newVal = (numChildren - 2.0)*dm.getDistance(child1->data.i,child2->data.i) - activeRowSums[child1->data.i] - activeRowSums[child2->data.i];
        if( newVal<minVal ){
          minVal = newVal;
          min_i = child1->data.i;
          min_j = child2->data.i;
        }
      }
    }

    //join the min pair into a parent node
    child1 = dm.getIdentifier(min_i);
    child2 = dm.getIdentifier(min_j);
    string_int si;
    si.i = vacantRows.back();
    vacantRows.pop_back();
    StringIntTree::Node *inner_parent = tree.detachFromParentAndAddAsSiblings(child1, child2, si);
    dm.setIdentifier(si.i,inner_parent);

    //UPDATE ONLY DISTANCES to inner_parent for the children of parent
    //(the distances to all other nodes will be computed later when the root is decided)
    //update the rowsums for the children of parent
    double r_2=(numChildren-2.0);
    double lambda = 0.5 + (varianceRowSums[min_j] - varianceRowSums[min_i])/(2.0*r_2*varianceMatrix.getDistance(min_i,min_j));
    if ( lambda > 1 ) lambda = 1.0;
    else if ( lambda < 0 ) lambda = 0;
    else if ( varianceMatrix.getDistance(min_i,min_j) == 0 )//i.e. lambda==NaN
      lambda = 0.5;
    
    double d1u = 0.5*(dm.getDistance(min_i,min_j) + (activeRowSums[min_i]-activeRowSums[min_j])/r_2);
    double d2u = 0.5*(dm.getDistance(min_i,min_j) + (activeRowSums[min_j]-activeRowSums[min_i])/r_2);

    StringIntTree::Node *child = parent->getRightMostChild();
    double rowSum = 0;
    double varianceSum = 0;
    for( ; child!=NULL ; child=child->getLeftSibling() ){
      if ( child==inner_parent ) continue;
      double dui = lambda*dm.getDistance(min_i,child->data.i) + (1.0-lambda)*dm.getDistance(min_j,child->data.i) - lambda*d1u - (1.0-lambda)*d2u;
      dm.setDistance(inner_parent->data.i,child->data.i, dui);
      double vui = lambda*varianceMatrix.getDistance(min_i,child->data.i) + (1-lambda)*varianceMatrix.getDistance(min_j,child->data.i)
        - lambda*(1-lambda)*varianceMatrix.getDistance(min_i,min_j);
      varianceMatrix.setDistance(inner_parent->data.i,child->data.i, vui);

      double dist2iandj = dm.getDistance(min_i,child->data.i) + dm.getDistance(min_j,child->data.i);
      double variance2iandj = varianceMatrix.getDistance(min_i,child->data.i) + varianceMatrix.getDistance(min_j,child->data.i);

      //update rowsums
      rowSum += dm.getDistance(inner_parent->data.i,child->data.i);
      varianceSum += varianceMatrix.getDistance(inner_parent->data.i,child->data.i);
      activeRowSums[child->data.i] += (dm.getDistance(inner_parent->data.i,child->data.i) -dist2iandj);
      varianceRowSums[child->data.i] += (varianceMatrix.getDistance(inner_parent->data.i,child->data.i) -variance2iandj);
    }
    activeRowSums[inner_parent->data.i] = rowSum;      
    varianceRowSums[inner_parent->data.i] = varianceSum;      
    dm.setDistance(inner_parent->data.i,inner_parent->data.i,0);
    varianceMatrix.setDistance(inner_parent->data.i,inner_parent->data.i,0);
    //----

  }//end inner NJ loop
  //----
}


StringIntTree
FNJ_collapse_method( String_DM &orig_dm ){


  numLeafs = orig_dm.getSize();
  //PRINT(orig_dm);

  //----
  //the tree that we are building
  //initiate it to a star.
  string_int strint;
  strint.i = -1;
  StringIntTree tree_tmp(strint);
  tree = tree_tmp;
  root = tree.getRoot();
  vector<StringIntTree::Node *> nodes(numLeafs);
  for( int i=0 ; i<numLeafs ; i++ ){
    string_int strint;
    strint.s = orig_dm.getIdentifier(i);
    strint.i = i;
    StringIntTree::Node *child = root->addChild(strint);
    nodes[i] = child;
  }
  //----


  //----
  //copy the distance matrix and compute the row sums
  dm.resize(2*numLeafs);
  varianceMatrix.resize(2*numLeafs);
  //the row sum vector contains the sum of distances between the
  //active nodes, i.e. the children of the root.
  activeRowSums.resize(2*numLeafs);
  varianceRowSums.resize(2*numLeafs);
  
  for( int i=0 ; i<numLeafs ; i++ ){
    dm.setIdentifier(i,nodes[i]);
    varianceMatrix.setIdentifier(i,nodes[i]);
    dm.setDistance(i,i,0);
    varianceMatrix.setDistance(i,i,0);
  }
  for( int i=numLeafs ; i<2*numLeafs ; i++ ){
    dm.setIdentifier(i,NULL);
  }

  for ( int i=0 ; i<numLeafs ; i++ ){
    double rowSum = 0;
    double varianceSum = 0;
    for( int j=0 ; j<numLeafs ; j++ ){
      if( j==i ) continue;
      double d = orig_dm.getDistance(i,j);
      rowSum += d;
      dm.setDistance(i,j,d);
      //d = d/(1.0*sequenceLength);
      varianceSum += d;
      varianceMatrix.setDistance(i,j,d); 
    } 
    
    activeRowSums[i] = rowSum;
    varianceRowSums[i] = varianceSum;
    //PRINT(i);PRINT(rowSum);
  }

  //the rows that are vacant
  vacantRows.clear();
  for( int i=2*numLeafs-1 ; i>numLeafs-1 ; i-- ){
    vacantRows.push_back(i);
  }
  //----

  //----
  //FNJ find the visible pairs
  // In the vector at position i we have the row id of the visible
  // sibling of the node in row i is
  vector<int> visibleSibling(2*numLeafs);
  for( int i=0 ; i<numLeafs ; i++ ){
    int minSibling = -1;
    double minVal =  FLT_MAX;
    for ( int j=0 ; j<numLeafs ; j++ ){
      if( j==i ) continue;
      double newVal = (numLeafs - 2.0)*dm.getDistance(i,j) - activeRowSums[i] - activeRowSums[j];
      if( newVal<minVal ){
        minVal = newVal;
        minSibling=j;
      }
    }

    visibleSibling[i] = minSibling;
  }
  
  //---
  //THE FNJ LOOP
  numChildrenOfRoot = numLeafs;
  for ( ; numChildrenOfRoot > 3 ; ) {
   
    //        LINE();tree.drawTree(cout);
    //find the min visible pair for the children of the root
    int min_i = -1;
    int min_j = -1;
    double minVal =  FLT_MAX;
    StringIntTree::Node *child = root->getRightMostChild();
    for( ; child!=NULL ; child=child->getLeftSibling() ){
      //PRINT(activeRowSums[child->data.i]);
      int visibleS = visibleSibling[child->data.i];
      double newVal = (numChildrenOfRoot - 2.0)*dm.getDistance(child->data.i,visibleS) - activeRowSums[child->data.i] - activeRowSums[visibleS];
      if( newVal<minVal ){
        minVal = newVal;
        min_i = child->data.i;
        min_j = visibleS;
      }
    }
    // PRINT(minVal);
    //PRINT(min_i);PRINT(min_j);


    //join the min pair into a parent node
    StringIntTree::Node *child1 = dm.getIdentifier(min_i);
    StringIntTree::Node *child2 = dm.getIdentifier(min_j);
    //cout << (void *) child1 <<endl;    cout << (void *) child2 <<endl;    
    int oldRootChild1 = child1->data.i;
    int oldRootChild2 = child2->data.i;
    string_int si;
    si.i = vacantRows.back();
    vacantRows.pop_back();
    StringIntTree::Node *parent = tree.detachFromParentAndAddAsSiblings(child1, child2, si);
    dm.setIdentifier(si.i,parent);
    numChildrenOfRoot--;
    
    //remove the distances to min_i and min_j from the activeRowSums of
    //the children of the root
    child = root->getRightMostChild();
    for( ; child!=NULL ; child=child->getLeftSibling() ){
      if( child==parent ) continue;
      double dist2iandj = dm.getDistance(min_i,child->data.i) + dm.getDistance(min_j,child->data.i);
      activeRowSums[child->data.i] -= dist2iandj;
      double variance2iandj = varianceMatrix.getDistance(min_i,child->data.i) + varianceMatrix.getDistance(min_j,child->data.i);
      varianceRowSums[child->data.i] -= variance2iandj;
    }
    //***************************
    //First NJ correction
    //Compute NJ on the 4 levels of parent and include the distance to
    //the other children of the root
    //collapse four levels
    if ( false ){
      collapse_children_levels(parent,COLLAPSE_LEVELS);
      if ( parent->getNumChildren() > 2 ){
        //compute the rowsums for the children of parent
        child1 = parent->getRightMostChild();
        for( ; child1!=NULL ; child1=child1->getLeftSibling() ){
          double rowSum = 0;
          double varianceSum = 0;
          //compute the sum to the distance to the children of parent
          child2 = parent->getRightMostChild();
          for( ; child2!=NULL ; child2=child2->getLeftSibling() ){
            if ( child2==child1 ) continue;
            rowSum += dm.getDistance(child1->data.i, child2->data.i);
            varianceSum += varianceMatrix.getDistance(child1->data.i, child2->data.i);
          }

          //compute the sum to the children of root
          child2 = root->getRightMostChild();
          for( ; child2!=NULL ; child2=child2->getLeftSibling() ){
            if ( child2==parent ) continue;
            rowSum += dm.getDistance(child1->data.i, child2->data.i);
            varianceSum += varianceMatrix.getDistance(child1->data.i, child2->data.i);
          }

          activeRowSums[child1->data.i] = rowSum;
          varianceRowSums[child1->data.i] = varianceSum;
        }

        //----
        run_regular_NJ_on_children(parent);
        //---
      }
    }
    //***********************************
    //Second NJ correction
    //Rebuild 4 levels without including the other children of the root.
    
    //collapse four levels 
    collapse_children_levels(parent,COLLAPSE_LEVELS);
    StringIntTree::NodeVector collapsedParentChildren;
    //        LINE();tree.drawSubtree(cout,parent);
    //compute the rowsums for the children of parent
    child1 = parent->getRightMostChild();
    for( ; child1!=NULL ; child1=child1->getLeftSibling() ){
      collapsedParentChildren.push_back(child1);
      double rowSum = 0;
      double varianceSum = 0;
      
      //compute the sum to the distance to the children of parent
      child2 = parent->getRightMostChild();
      for( ; child2!=NULL ; child2=child2->getLeftSibling() ){
        if ( child2==child1 ) continue;
        rowSum += dm.getDistance(child1->data.i, child2->data.i);
        varianceSum += varianceMatrix.getDistance(child1->data.i, child2->data.i);
      }

      activeRowSums[child1->data.i] = rowSum;
      varianceRowSums[child1->data.i] = varianceSum;
    }
    if( parent->getNumChildren()>3 ){

      //REBUILD using NJ without using rest of root children
      if( SUBTREE_METHOD=='W' )
        run_subtree_WeighborExternally_on_children(parent);
      else
        run_subtree_NJ_on_children(parent);
      assert(parent->getNumChildren()<=3);
    }
    //********************
    //FIND root using NJ and recompute the distances from the new
    //nodes to all other nodes.  This time all distances need to
    //included.
    //compute the rowsums for the collapsed children of parent
    if( parent->getNumChildren()>2 ){
      for( size_t i=0 ; i<collapsedParentChildren.size() ; i++ ){
        child1 = collapsedParentChildren[i];
        double rowSum = 0;
        double varianceSum = 0;
        
        //compute the sum of distances to the collapsed children of parent
        for( size_t j=0 ; j<collapsedParentChildren.size() ; j++ ){
          child2 = collapsedParentChildren[j];
          if ( child2==child1 ) continue;
          rowSum += dm.getDistance(child1->data.i, child2->data.i);
          varianceSum += varianceMatrix.getDistance(child1->data.i, child2->data.i);
        }

        //compute the sum of distances to the children of root
        child2 = root->getRightMostChild();
        for( ; child2!=NULL ; child2=child2->getLeftSibling() ){
          if ( child2==parent ) continue;
          rowSum += dm.getDistance(child1->data.i, child2->data.i);
          varianceSum += varianceMatrix.getDistance(child1->data.i, child2->data.i);
        }

        activeRowSums[child1->data.i] = rowSum;
        varianceRowSums[child1->data.i] = varianceSum;
      }

      //---
      StringIntTree parent_tree;
      tree.splittTree(parent,parent_tree);//takes linear time
      //      LINE();parent_tree.drawTree(cout);
      //start NJ to find the root, i.e. join nodes until no more nodes
      //can be joined. Infere a node on the last edge and connect it to
      //the root of the original tree.
      //SEPARATOR();parent_tree.drawTree(cout);
      while( collapsedParentChildren.size()>2 ){
        //        SEPERATOR();
        //find the  minimal sibling pair
        minVal = FLT_MAX;
        StringIntTree::Node *min_parent = NULL;
        size_t collapsed_index_i = 1000000;
        size_t collapsed_index_j = 1000000;

        for( size_t i=0 ; i<collapsedParentChildren.size() ; i++ ){
          child1 = collapsedParentChildren[i];//LINE();PRINT("checking");PRINT(child1->getNodeId());
          StringIntTree::Node *inner_parent = NULL;
          //find sibling
          for( size_t j=i+1 ; j<collapsedParentChildren.size() ; j++ ){
            child2 = collapsedParentChildren[j];
            //LINE();PRINT(child1->getNodeId());PRINT(child2->getNodeId());PRINT(j);PRINT(collapsedParentChildren.size());
            if( child1==child2 ) continue;
            StringIntTree::Node *c1p = child1->getParent();
            StringIntTree::Node *c2p = child2->getParent();
            if(  c1p==c2p )
              inner_parent = c1p;
            else if ( c1p!=NULL && c1p->getParent() == child2 )
              inner_parent = c1p;
            else if( c2p!=NULL && c2p->getParent() == child1 )
              inner_parent = c2p;
            else //not siblings
              continue;
            //  LINE();PRINT("CHECKING");PRINT(child1->getNodeId());PRINT(child2->getNodeId());
            //child1 and child2 are siblings and inner_parent the parent
            double newVal = (collapsedParentChildren.size()+(numChildrenOfRoot-1) - 2.0)*dm.getDistance(child1->data.i,child2->data.i)
              - activeRowSums[child1->data.i] - activeRowSums[child2->data.i];
            //PRINT(newVal);PRINT(activeRowSums[child1->data.i]);PRINT(activeRowSums[child2->data.i]);            
            if( newVal<minVal ){
              minVal = newVal;
              min_i = child1->data.i;
              min_j = child2->data.i;
              //LINE();PRINT(child1->getNodeId());PRINT(child2->getNodeId());
              collapsed_index_i = i;
              collapsed_index_j = j;
              //   min_i = child1->data.i;
              //               min_j = child2->data.i;
              //               //LINE();PRINT(child1->getNodeId());PRINT(child2->getNodeId());
              //               collapsed_index_i = i;
              //               collapsed_index_j = j;
              min_parent = inner_parent;
            }
          }
        }
        
        //insert min_parent into the collapsed vector and remove its children
        if( collapsed_index_j!=collapsedParentChildren.size()-1 ){
          if ( collapsed_index_i!=collapsedParentChildren.size()-1 )
            collapsedParentChildren[collapsed_index_j] = collapsedParentChildren.back();
          else
            collapsed_index_i = collapsed_index_j;
        }     
        collapsedParentChildren.pop_back();
        collapsedParentChildren[collapsed_index_i] = min_parent;
        // LINE();
        //         PRINT(min_i);
        //         PRINT((void*)dm.getIdentifier(min_i));
        //         PRINT(dm.getIdentifier(min_i)->getNodeId());PRINT(dm.getIdentifier(min_j)->getNodeId());PRINT(min_parent->getNodeId());
        //        LINE();
        //----
        //update distances to min_parent
        double r_2 = (collapsedParentChildren.size()+1+(numChildrenOfRoot-1) - 2.0);
        double lambda = 0.5 + (varianceRowSums[min_j] - varianceRowSums[min_i])/(2.0*r_2
                                                                                 *varianceMatrix.getDistance(min_i,min_j));
        if ( lambda > 1 ) lambda = 1.0;
        else if ( lambda < 0 ) lambda = 0;
        else if ( varianceMatrix.getDistance(min_i,min_j) == 0 )//i.e. lambda==NaN
          lambda = 0.5;
    
        double d1u = 0.5*(dm.getDistance(min_i,min_j) + (activeRowSums[min_i]-activeRowSums[min_j])/r_2);
        double d2u = 0.5*(dm.getDistance(min_i,min_j) + (activeRowSums[min_j]-activeRowSums[min_i])/r_2);
      
        for( int i=0 ; i<2*numLeafs ; i++ ){
          if ( dm.getIdentifier(i) == NULL ) continue;
          double dui = lambda*dm.getDistance(min_i,i) + (1.0-lambda)*dm.getDistance(min_j,i) - lambda*d1u - (1.0-lambda)*d2u;
          dm.setDistance(min_parent->data.i,i, dui);
          double vui = lambda*varianceMatrix.getDistance(min_i,i) + (1-lambda)*varianceMatrix.getDistance(min_j,i) - lambda*(1-lambda)*varianceMatrix.getDistance(min_i,min_j);
          varianceMatrix.setDistance(min_parent->data.i,i, vui);
        }
        dm.setDistance(min_parent->data.i,min_parent->data.i,0);
      
        //compute the row sums of min parent 
        double parentRowSum =0;
        double parentVarianceSum = 0;
        child = root->getRightMostChild();
        for( ; child!=NULL ; child=child->getLeftSibling() ){
          parentRowSum += dm.getDistance(min_parent->data.i,child->data.i);
          parentVarianceSum += varianceMatrix.getDistance(min_parent->data.i,child->data.i);
        }
        for( size_t i=0 ; i<collapsedParentChildren.size() ; i++ ){
          child = collapsedParentChildren[i];
          if( child==min_parent ) continue;

          double dist2iandj = dm.getDistance(min_i,child->data.i) + dm.getDistance(min_j,child->data.i);
          double variance2iandj = varianceMatrix.getDistance(min_i,child->data.i) + varianceMatrix.getDistance(min_j,child->data.i);
          //PRINT(dist2iandj);
          parentRowSum += dm.getDistance(min_parent->data.i,child->data.i);
          parentVarianceSum += varianceMatrix.getDistance(min_parent->data.i,child->data.i);
          activeRowSums[child->data.i] += (dm.getDistance(min_parent->data.i,child->data.i) - dist2iandj);
          varianceRowSums[child->data.i] += (varianceMatrix.getDistance(min_parent->data.i,child->data.i) - variance2iandj);
          
          //PRINT(activeRowSums[child->data.i]);PRINT(child->getNodeId());
        }    
        activeRowSums[min_parent->data.i] = parentRowSum;
        varianceRowSums[min_parent->data.i] = parentVarianceSum;
        
        //PRINT(parentRowSum);PRINT(min_parent->data.i);
      }
      //      LINE();parent_tree.drawTree(cout);
      ASSERT_EQ(collapsedParentChildren.size(),2);
      child1 = collapsedParentChildren[0];
      child2 = collapsedParentChildren[1];
      //    LINE();PRINT(child1->getNodeId());PRINT(child2->getNodeId());
      assert(child1->getParent()==child2 || child2->getParent()==child1);
      si.i = vacantRows.back();
      vacantRows.pop_back();
      if( child1->getParent()==child2 )
        parent = parent_tree.insertNodeOnPathToParent(child1, si);
      else if( child2->getParent()==child1 )
        parent = parent_tree.insertNodeOnPathToParent(child2, si);
      else
        PROG_ERROR("Shouldn't come here... the last two nodes aren't neighbors");
      dm.setIdentifier(si.i, parent);

      parent_tree.reRootAt(parent);
      //      LINE();parent_tree.drawTree(cout);
      tree.joinTreeAtRoot(parent_tree);
    }
    //*******************
      
    
    //UPDATE ALL DISTANCES to parent
    assert( parent->getNumChildren() == 2 );
    min_i = parent->getRightMostChild()->data.i;
    min_j = parent->getRightMostChild()->getLeftSibling()->data.i;
    //lambda the BIONJ way
    double lambda = 0.5 + (varianceRowSums[min_j] - varianceRowSums[min_i])/(2.0*(numChildrenOfRoot+1-2.0)*varianceMatrix.getDistance(min_i,min_j));
    //PRINT(lambda);PRINT(varianceMatrix.getDistance(min_i,min_j));PRINT(dm.getDistance(min_i,min_j));
    if ( lambda > 1 ) lambda = 1.0;
    else if ( lambda < 0 ) lambda = 0;
    else if ( varianceMatrix.getDistance(min_i,min_j) == 0 )//i.e. lambda==NaN
      lambda = 0.5;
    
    double d1u = 0.5*(dm.getDistance(min_i,min_j) + (activeRowSums[min_i]-activeRowSums[min_j])/(numChildrenOfRoot+1-2.0));
    double d2u = 0.5*(dm.getDistance(min_i,min_j) + (activeRowSums[min_j]-activeRowSums[min_i])/(numChildrenOfRoot+1-2.0));
      
    for( int i=0 ; i<2*numLeafs ; i++ ){
      if ( dm.getIdentifier(i) == NULL ) continue;
      double dui = lambda*dm.getDistance(min_i,i) + (1.0-lambda)*dm.getDistance(min_j,i) - lambda*d1u - (1.0-lambda)*d2u;
      dm.setDistance(parent->data.i,i, dui);
      double vui = lambda*varianceMatrix.getDistance(min_i,i) + (1-lambda)*varianceMatrix.getDistance(min_j,i) - lambda*(1-lambda)*varianceMatrix.getDistance(min_i,min_j);
      varianceMatrix.setDistance(parent->data.i,i, vui);
    }
    dm.setDistance(parent->data.i,parent->data.i,0);
    varianceMatrix.setDistance(parent->data.i,parent->data.i,0);

    //update the row sums
    double parentRowSum =0;
    double parentVarianceSum =0;
    child = root->getRightMostChild();
    for( ; child!=NULL ; child=child->getLeftSibling() ){
      if ( child==parent ) continue;
      double dist2iandj = dm.getDistance(min_i,child->data.i) + dm.getDistance(min_j,child->data.i);
      double variance2iandj = varianceMatrix.getDistance(min_i,child->data.i) + varianceMatrix.getDistance(min_j,child->data.i);
      //PRINT(dist2iandj);
      parentRowSum += dm.getDistance(parent->data.i,child->data.i);
      parentVarianceSum += varianceMatrix.getDistance(parent->data.i,child->data.i);
      activeRowSums[child->data.i] += (dm.getDistance(parent->data.i,child->data.i));//the distance to the old min_i and min_j has already been removed
      varianceRowSums[child->data.i] += (varianceMatrix.getDistance(parent->data.i,child->data.i));//the distance to the old min_i and min_j has already been removed
      // PRINT(activeRowSums[child->data.i]);
    }
    
    activeRowSums[parent->data.i] = parentRowSum;
    varianceRowSums[parent->data.i] = parentVarianceSum;

    //FIX VISIBLE PAIRS
    //go through and find all vissible pairs containing the old two nodes
    //and exchange the old nodes for the new parent node
    child = root->getRightMostChild();
    ASSERT_EQ(root->getNumChildren(),numChildrenOfRoot);
    //find the visible pair of parent
    int min_parent_sib=-1;
    minVal =  FLT_MAX;
    for( ; child!=NULL ; child=child->getLeftSibling() ){
      if( child==parent) continue;
      double newVal = (numChildrenOfRoot - 2.0)*dm.getDistance(child->data.i,parent->data.i) - activeRowSums[child->data.i] - activeRowSums[parent->data.i];
      if( newVal<minVal ){
        minVal = newVal;
        min_parent_sib = child->data.i; 
      }

      int visibleS = visibleSibling[child->data.i];
      if( visibleS==oldRootChild1 || visibleS==oldRootChild2 )
        visibleSibling[child->data.i] = parent->data.i;
    }
    visibleSibling[parent->data.i] = min_parent_sib;
    //        PRINT(oldRootChild1);PRINT(oldRootChild2);PRINT(parent->data.i);PRINT(min_parent_sib);
  }//end of FNJ loop
  //---------


  //-----------
  //COLLAPSE THE ROOT four levels and recompute
  //collapse four levels 
  collapse_children_levels(root,COLLAPSE_LEVELS);
  //compute the rowsums for the children of parent
  StringIntTree::Node *child1 = root->getRightMostChild();
  for( ; child1!=NULL ; child1=child1->getLeftSibling() ){
    double rowSum = 0;

    //compute the sum to the distance to the children of root
    StringIntTree::Node *child2 = root->getRightMostChild();
    for( ; child2!=NULL ; child2=child2->getLeftSibling() ){
      if ( child2==child1 ) continue;
      rowSum += dm.getDistance(child1->data.i, child2->data.i);
    }
    activeRowSums[child1->data.i] = rowSum;
  }

  //----
  if( SUBTREE_METHOD=='W' )
    run_subtree_WeighborExternally_on_children(root);
  else
    run_subtree_NJ_on_children(root);

  //---
  //  LINE();tree.drawTree(cout);
  return tree;
}




void
print_usage(){
  cout << "usage: FNJ_collapse_method -i dmatrix -m 20"  <<endl;
  cout << "-s  subtree method options \'B\'=BioNJ \'W\'=Weighbor, default B"<<endl;
  cout << "-seqlen sequence length (only needed if Weighbor is used"<<endl;
  cout << "-cl number of levels to collaps, default 3."<<endl;
  exit(1);
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

  sequenceLength = -1;
  char *seqlen = GET_OPTION_VAL("-seqlen");
  if ( seqlen != NULL)
    sequenceLength = atoi(seqlen);


  SUBTREE_METHOD = 'B';
  char *stmstr = GET_OPTION_VAL("-s");
  if( stmstr!=NULL ){
    if( stmstr[0]=='W' ){
      SUBTREE_METHOD='W';
      if( sequenceLength<=0 )
        print_usage();
    }
    print_usage();
  }

  COLLAPSE_LEVELS=3;
  char *cstr = GET_OPTION_VAL("-cl");
  if( cstr!=NULL ){
    COLLAPSE_LEVELS = atoi(cstr);
  }

  StringIntTree::NodeVector nodes;
  nodes.reserve(400);
  for ( int i=0 ; i<m ; i++ ){
    DM dm(in);

    StringIntTree tree = FNJ_collapse_method(dm);

    tree.addNodesInPostfixOrder(nodes);
    for( size_t j=0 ; j<nodes.size() ; j++ )
      nodes[j]->data.i=-1;
    nodes.clear();
    cout << tree << endl;
  }
  in.close();

return 1;
}
