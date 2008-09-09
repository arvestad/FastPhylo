//--------------------------------------------------
//                                        
// File: NeighborJoining.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: NeighborJoining.hpp,v 1.28 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef NEIGHBORJOINING_HPP
#define NEIGHBORJOINING_HPP

#include "DistanceMatrix.hpp"
#include <string>
#include <math.h>
#include "log_utils.hpp"
#include <iostream>
#include <float.h>
#include "SequenceTree.hpp"
#include <algorithm>

//-------------------------- NEIGHBOR METHODS --------------------------------
//
// There are three implemented methods: the regular NJ algorithm by
// Saitou and Nei, the Fast NJ algorithm by Elias and Lagergren, and
// the Bio NJ algorithm by Gascuel.
//
// There are two approaches to reconstructing a tree, either to use
// the general method with a distance method as input or to use the
// specific methods in which the identifiers are nodes connected in a
// star. The specific methods thereby allows for resolving a star
// which is a subtree of any tree.
//
enum NJ_method{
  NJ, 
  FNJ,
  BIONJ
};


//-------------------------------------------------------------------------
// GENERAL METHOD
// A tree is built from the distance matrix using the supplied method.
// 
//
void
computeNJTree(StrDblMatrix &dm, SequenceTree &resultTree, NJ_method m=NJ );

//---------------------- NEIGHBOR JOINING ----------------------------------------
//Takes a distance matrix in which the identifiers are tree nodes.
//The tree nodes should be in connected in a star. That may be a subtree
//in a big tree.

template <class TreeNode_type, class Data>
void
computeNeighborJoiningTree( DistanceMatrix< TreeNode_type *, double, 
			    Data_init<TreeNode_type*>, Data_printOn<TreeNode_type *>, 
			    Data_init<double>, Data_printOn<double> > &dm, 
			    Data defaultNodeData ){

  const size_t origNumNodes = dm.getSize(); 
#ifndef NDEBUG
  //make sure that the input nodes are in a star  
  TreeNode_type *starroot = dm.getIdentifier(0)->getParent();
  for ( size_t i = 0 ; i < origNumNodes ; i++ )
    if ( dm.getIdentifier(i)->getParent() != starroot ){
      PROG_ERROR("The input nodes to NJ are not connected in a star");
    }
#endif
  //compute the row sums
  double rowSums[origNumNodes];
  for ( size_t row = 0 ; row < origNumNodes ; row++ ){
    double sum = 0;
    size_t i = 0;
    for (  ; i < origNumNodes ; i++ ){
      double d=dm.getDistance(row,i);
      if(isfinite(d))
	sum += d;
      else{
	USER_ERROR("Distance Matrix contains a non finite number: " << d);
      }
    }
    rowSums[row] = sum;
  }


  //----------------
  // NJ ITERATION
  size_t numNodes = dm.getSize();
  while ( numNodes > 3 ) {
    assert(dm.getSize() == numNodes);
    //find the minimal value
    double minVal = FLT_MAX;
    size_t mini = 1000000;
    size_t minj = 1000000;
    for ( size_t i = 0 ; i < numNodes ; i++ ){
      for ( size_t j = i+1 ; j < numNodes ; j++ ){
        double newVal = (numNodes - 2.0)*dm.getDistance(i,j) - rowSums[i] - rowSums[j];
        if ( newVal < minVal ){
          minVal = newVal;
          mini = i;
          minj = j;
        }
      }
    }

    //make sure that minj is the last row in the matrix
    if ( minj != numNodes -1 ){
      if(mini==numNodes-1){
        mini=minj;
        minj=numNodes-1;
        double tmp = rowSums[numNodes-1];
        rowSums[numNodes-1] = rowSums[mini];
        rowSums[mini] = tmp;
      }else{
	//        PRINT(minj);PRINT(numNodes);
        dm.swapRowToLast(minj);
	std::swap(rowSums[numNodes-1],rowSums[minj]);
        minj = numNodes-1;
      }
    }

    TreeNode_type *newparent = dm.getIdentifier(mini)->getTree()->
      detachFromParentAndAddAsSiblings(dm.getIdentifier(mini),dm.getIdentifier(minj), defaultNodeData);
    dm.setIdentifier(mini, newparent);
    
    // UPDATE DISTANCES
    for ( size_t i = 0 ; i < numNodes-1 ; i++ ){//skip last row
      double dist2iandj = dm.getDistance(mini,i) + dm.getDistance(minj,i);
      // regular nj update function:
      dm.setDistance(mini,i, dist2iandj * 0.5); 
    
      //update rowsums
      rowSums[i] = rowSums[i] - dist2iandj + dm.getDistance(mini,i);
    }

    //remove the last row of the matrix
    dm.removeLastRow();
    numNodes--;    

    //recompute the row sum for the parent
    dm.setDistance(mini,mini,0);
    double sum = 0;
    for ( size_t i = 0 ; i < numNodes ; i++ )
      sum += dm.getDistance(mini,i);
    rowSums[mini] = sum;

  }
  // END ITERATION
  //--------------
}


//-------------------- BIO NJ -------------------------------------------------------------
//
//

template <class TreeNode_type, class Data>
void
computeBioNJTree( DistanceMatrix< TreeNode_type *, double, 
		  Data_init<TreeNode_type*>, Data_printOn<TreeNode_type *>, 
		  Data_init<double>, Data_printOn<double> > &dm, 
		  Data defaultNodeData ){
  
  const size_t origNumNodes = dm.getSize();
  //make sure that the input nodes are in a star  
#ifndef NDEBUG
  TreeNode_type *starroot = dm.getIdentifier(0)->getParent();
  for ( size_t i = 0 ; i < origNumNodes ; i++ ){
    if ( dm.getIdentifier(i)->getParent() != starroot ){
      PROG_ERROR("The input nodes to NJ are not connected in a star");
    }
  }
#endif

  //the variance matrix
  StrDblMatrix V(origNumNodes);//the variance is copied from the regular distance matrix
  for(size_t i=0;i<origNumNodes;i++){
    Sequence_double data;
    for(size_t j=i;j<origNumNodes;j++){
      V.setDistance(i,j,dm.getDistance(i,j));
    }
  }

  //compute the row sums
  double rowSums[origNumNodes];
  double varianceRowSums[origNumNodes];
  for ( size_t row = 0 ; row < origNumNodes ; row++ ){
    double sum = 0;
    double sumV =0;
    size_t i = 0;
    for (  ; i < origNumNodes ; i++ ){
      double d=dm.getDistance(row,i);
      if(isfinite(d)){
	sum += d;
	sumV += V.getDistance(row,i);
      }else{
	USER_ERROR("Distance Matrix contains a non finite number: " << d);
      }
    }
    rowSums[row] = sum;
    varianceRowSums[row] = sum;
  }


  //----------------
  //ITERATION
  size_t numNodes = dm.getSize();
  while (  numNodes > 3) {
    assert(dm.getSize() == numNodes);
    //find the minimal value
    double minVal = FLT_MAX;
    size_t mini = 9999999;
    size_t minj = 9999999;
    for ( size_t i = 0 ; i < numNodes ; i++ ){
      for ( size_t j = i+1 ; j < numNodes ; j++ ){
        double newVal = (numNodes - 2.0)*dm.getDistance(i,j) - rowSums[i] - rowSums[j];
        if ( newVal < minVal ){
          minVal = newVal;
          mini = i;
          minj = j;
        }
      }
    }

    //make sure that minj is the last row in the matrix
    if ( minj != numNodes -1 ){
      if(mini==numNodes-1){
        mini=minj;
        minj=numNodes-1;
        double tmp = rowSums[numNodes-1];
        rowSums[numNodes-1] = rowSums[mini];
        rowSums[mini] = tmp;
	tmp = varianceRowSums[numNodes-1];
        varianceRowSums[numNodes-1] = varianceRowSums[mini];
        varianceRowSums[mini] = tmp;
      }else{
	//        PRINT(minj);PRINT(numNodes);
        dm.swapRowToLast(minj);
	V.swapRowToLast(minj);
	std::swap(rowSums[numNodes-1],rowSums[minj]);
	std::swap(varianceRowSums[numNodes-1],varianceRowSums[minj]);
        minj = numNodes-1;
      }
    }


    //COMPUTE lamba
    double V_a2b = V.getDistance(mini,minj);
    double lambda = 0.5 + (varianceRowSums[minj] - varianceRowSums[mini])/(2*(V.getSize()-2)*V_a2b);
    if ( lambda > 1 ) lambda = 1.0;
    else if ( lambda < 0 ) lambda = 0;
    else if ( V_a2b == 0 )//i.e. lambda==NaN
      lambda = 0.5;
    
    
   double D_a2b = dm.getDistance(mini,minj);
    
   TreeNode_type *a = dm.getIdentifier(mini);
   TreeNode_type *b = dm.getIdentifier(minj);
   
   TreeNode_type *newparent = dm.getIdentifier(mini)->getTree()->
     detachFromParentAndAddAsSiblings(dm.getIdentifier(mini),dm.getIdentifier(minj), defaultNodeData);
    dm.setIdentifier(mini, newparent);
    
    double D_a2parent = 0.5*(D_a2b+(rowSums[mini]-rowSums[minj])/(numNodes-2));
    double D_b2parent = 0.5*(D_a2b+(rowSums[minj]-rowSums[mini])/(numNodes-2));
    EDGE(a)=D_a2parent;
    EDGE(b)=D_b2parent;
  

    // UPDATE DISTANCES
    for ( size_t i = 0 ; i < numNodes-1 ; i++ ){//skip last row
      double dist2ab = dm.getDistance(mini,i) + dm.getDistance(minj,i);
      double newdist = lambda*dm.getDistance(mini,i)+(1-lambda)*dm.getDistance(minj,i) 
	- lambda*D_a2parent - (1-lambda)*D_b2parent;
      
      dm.setDistance(mini,i,newdist);
      rowSums[i] = rowSums[i] - dist2ab + dm.getDistance(mini,i);

      //Update variance
      double var2ab = V.getDistance(i,mini)+V.getDistance(i,minj);
      V.setDistance(i,mini, lambda*V.getDistance(i,mini) + (1-lambda)*V.getDistance(i,minj) - lambda*(1-lambda)*V_a2b);
      varianceRowSums[i] = varianceRowSums[i] - var2ab + V.getDistance(i,mini);
    }

    //remove the last row of the matrix
    dm.removeLastRow();
    V.removeLastRow();
    numNodes--;
    //recompute the row sum for the parent
    dm.setDistance(mini,mini,0);
    V.setDistance(mini,mini,0);
    double sum = 0;
    double sumV = 0;
    for ( size_t i = 0 ; i < numNodes ; i++ ){
      sum += dm.getDistance(mini,i);
      sumV += V.getDistance(mini,i);
    }
    rowSums[mini] = sum;
    varianceRowSums[mini] = sumV;
    

  }
  // END ITERATION
  //--------------

  EDGE(dm.getIdentifier(0)) = dm.getDistance(0,1) + dm.getDistance(0,2)
    - 2 * dm.getDistance(1,2);
  EDGE(dm.getIdentifier(1)) = dm.getDistance(0,1) + dm.getDistance(1,2)
    - 2 * dm.getDistance(0,2);
  EDGE(dm.getIdentifier(2)) = dm.getDistance(0,2) + dm.getDistance(1,2)
    - 2 * dm.getDistance(0,1);
}


//------------------------------- FAST NEIGHBOR JOINING -----------------------
//
//
//



template <class TreeNode_type, class Data>
void
computeFNJTree( DistanceMatrix< TreeNode_type *, double, 
		Data_init<TreeNode_type*>, Data_printOn<TreeNode_type *>, 
		Data_init<double>, Data_printOn<double> > &dm, Data defaultNodeData ){

  const size_t origNumNodes = dm.getSize(); 

  //make sure that the input nodes are in a star  
#ifndef NDEBUG
  TreeNode_type *starroot = dm.getIdentifier(0)->getParent();
  for ( size_t i = 0 ; i < origNumNodes ; i++ )
    if ( dm.getIdentifier(i)->getParent() != starroot ){
      PROG_ERROR("The input nodes to NJ are not connected in a star");
    }
#endif

  //compute the row sums
  double rowSums[origNumNodes];
  for ( size_t row = 0 ; row < origNumNodes ; row++ ){
    double sum = 0;
    size_t i = 0;
    for (  ; i < origNumNodes ; i++ ){
      double d=dm.getDistance(row,i);
      if(isfinite(d))
	sum += d;
      else{
	USER_ERROR("Distance Matrix contains a non finite number: " << d);
      }
    }
    rowSums[row] = sum;
  }
  //compute visible set.
  size_t visible_set[origNumNodes];
  for ( size_t i = 0 ; i < origNumNodes ; i++ ){
    double minVal = FLT_MAX;
    size_t minNeigh = i;
    for ( size_t j = 0 ; j < i ; j++ ){
      double newVal = (origNumNodes - 2.0)*dm.getDistance(i,j) - rowSums[i] - rowSums[j];
      if ( newVal < minVal ){
	minVal = newVal;
	minNeigh = j;
      }
    }
    
    for ( size_t j = i+1 ; j < origNumNodes ; j++ ){
      double newVal = (origNumNodes - 2.0)*dm.getDistance(i,j) - rowSums[i] - rowSums[j];
      if ( newVal < minVal ){
	minVal = newVal;
	minNeigh = j;
      }
    }
    visible_set[i] = minNeigh;
    assert(minNeigh != i );
  }


  //----------------
  // NJ ITERATION
  size_t numNodes = dm.getSize();
  while ( numNodes > 3 ) {
    assert(dm.getSize() == numNodes);
    //find the minimal value
    size_t mini = 0;
    size_t minj = visible_set[0];
    double minVal = (numNodes - 2.0)*dm.getDistance(0,minj) - rowSums[0] - rowSums[minj];
    // min over visible set
    for ( size_t i = 1 ; i < numNodes ; i++ ){
      size_t j = visible_set[i];
      double newVal = (numNodes - 2.0)*dm.getDistance(i,j) - rowSums[i] - rowSums[j];
      if ( newVal < minVal ){
	minVal = newVal;
	mini = i;
	minj = j;
      }
    }

   
    //make sure that minj is the last row in the matrix
    if ( minj != numNodes -1 ){
      if(mini==numNodes-1){
        mini=minj;
        minj=numNodes-1;
        double tmp = rowSums[numNodes-1];
        rowSums[numNodes-1] = rowSums[mini];
        rowSums[mini] = tmp;

	//update the visible set so that no one points to numNodes-1
	for(size_t vi = 0; vi<numNodes; vi++)
	  if(visible_set[vi]==minj)
	    visible_set[vi] = mini;//those that point to 
      }
      else{// minj needs to be swapped to the last row

	//update visible set
	for(size_t vi = 0; vi<numNodes; vi++)
	  if(visible_set[vi]==numNodes-1)
	    visible_set[vi] = minj;//those that point to the last row before the swap
	  else if (visible_set[vi]==minj)
	    visible_set[vi] = mini;//those that point minj before the swap
	  
        dm.swapRowToLast(minj);
	std::swap(rowSums[numNodes-1],rowSums[minj]);
	std::swap(visible_set[numNodes-1],visible_set[minj]);
	minj = numNodes-1;

      }
    }
    else{//if minj == numNodes -1 
      for(size_t vi = 0; vi<numNodes; vi++)
	if(visible_set[vi]==minj)
	  visible_set[vi] = mini;//those that point to minj
    }
    
    //join the tree nodes    
    TreeNode_type *newparent = dm.getIdentifier(mini)->getTree()->
      detachFromParentAndAddAsSiblings(dm.getIdentifier(mini),dm.getIdentifier(minj), defaultNodeData);
    dm.setIdentifier(mini, newparent);


    // UPDATE DISTANCES
    for ( size_t i = 0 ; i < numNodes-1 ; i++ ){//skip last row
      double dist2iandj = dm.getDistance(mini,i) + dm.getDistance(minj,i);
      // regular nj update function:
      dm.setDistance(mini,i, dist2iandj * 0.5); 
    
      //update rowsums
      rowSums[i] = rowSums[i] - dist2iandj + dm.getDistance(mini,i);
    }

    //remove the last row of the matrix
    dm.removeLastRow();numNodes--;
    //*** update visible set.
    //recompute the row sum for the parent
    dm.setDistance(mini,mini,0);
    double sum = 0;
    for ( size_t i = 0 ; i < numNodes ; i++ )
      sum += dm.getDistance(mini,i);
    rowSums[mini] = sum;
    //compute visible neighbor
    minVal = FLT_MAX;
    size_t minNeigh = mini;
    for ( size_t j = 0 ; j < mini ; j++ ){
      double newVal = (numNodes - 2.0)*dm.getDistance(mini,j) - rowSums[mini] - rowSums[j];
      if ( newVal < minVal ){
	minVal = newVal;
	minNeigh = j;
      }
    }
    for ( size_t j = mini+1 ; j < numNodes ; j++ ){
      double newVal = (numNodes - 2.0)*dm.getDistance(mini,j) - rowSums[mini] - rowSums[j];
      if ( newVal < minVal ){
	minVal = newVal;
	minNeigh = j;
      }
    }
    
    visible_set[mini] = minNeigh;
  }
  // END ITERATION
  //--------------
}



#endif // NEIGHBORJOINING_HPP














