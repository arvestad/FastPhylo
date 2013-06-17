#include <float.h>
#include <math.h>
#include "NeighborJoining.hpp"
#include "Sequence.hpp"

typedef DistanceMatrix<SequenceTree::Node *,double,Data_init<SequenceTree::Node *>,Data_printOn<SequenceTree::Node *>,Data_init<double>,Data_printOn<double> > NJMatrix;
typedef FloatDistanceMatrix<SequenceTree::Node *,float,Data_init<SequenceTree::Node *>,Data_printOn<SequenceTree::Node *>,Data_init<float>,Data_printOn<float> > NJFloMatrix;


void
computeNJTree(StrDblMatrix &dm, SequenceTree &tree, NJ_method m ){

  //create a star tree
  Sequence_double defdata;
  defdata.dbl=-1;
  tree = SequenceTree(defdata);
  
  NJMatrix njdm(dm.getSize());
  
  for(size_t i=0;i<dm.getSize();i++){
    Sequence_double data;
    data.dbl = -1;
    data.s.name = dm.getIdentifier(i);
    SequenceTree::Node *node = tree.getRoot()->addChild(data);
    njdm.setIdentifier(i,node);
    for(size_t j=i;j<dm.getSize();j++){
      njdm.setDistance(i,j,dm.getDistance(i,j));
    }
  }

  //create NJ tree
  switch(m){
  case NJ:computeNeighborJoiningTree(njdm,defdata); break;
  case BIONJ: computeBioNJTree(njdm,defdata); break;
  case FNJ: computeFNJTree(njdm,defdata); break;
  default:
    PROG_ERROR("Unexpected method");
  }
}

// Mehmood's addition here
void
computeNJTree(StrFloMatrix &dm, SequenceTree &tree, NJ_method m ){
  //Rewrote this method because it unneccessarily created two distance matrices
  int rows = dm.getRows();
  int columns = dm.getColumns();

  int helpMatrixSize = 0;
  //Because DistanceMatrix copies automatically the distance (i,j) and (j,i) and so I used a vector of Rows
  std::vector<StrFloRow> helpMatrix;
  helpMatrix.resize(helpMatrixSize+1);
  helpMatrix[helpMatrixSize].resize(helpMatrixSize+1);

  for(int i = rows-1; i >= 0; --i){
  	++helpMatrixSize;
  	helpMatrix[helpMatrixSize-1].setIdentifier(dm.getIdentifier(i));
  	for(int j = i; j < columns; ++j){
  		helpMatrix[helpMatrixSize-1].setDistance(j-i, dm.getDistance(i,j));
  	}
  	dm.removeLastRow();
  	helpMatrix.resize(helpMatrixSize+1);
  	helpMatrix[helpMatrixSize].resize(helpMatrixSize+1);

  }

  //create a star tree
  Sequence_double defdata;
  defdata.dbl=-1;
  tree = SequenceTree(defdata);

  NJFloMatrix njflodm(1, columns);
  int rowCount = 0;
  for(int i=rows-1;i>=0;--i){
    Sequence_double data;
    data.dbl = -1;
    data.s.name = helpMatrix[i].getIdentifier();
    SequenceTree::Node *node = tree.getRoot()->addChild(data);
    njflodm.setIdentifier(rowCount,node);

    for(size_t j=rowCount;j<columns;j++){
    	float test = helpMatrix[i].getDistance(j-rowCount);
      	njflodm.setDistance(rowCount,j,test);
    }
  	--helpMatrixSize;
  	helpMatrix.resize(helpMatrixSize);
    njflodm.addRow();

    ++rowCount;
  }

  //create a star tree
/*  Sequence_double defdata;
  defdata.dbl=-1;
  tree = SequenceTree(defdata);

  NJMatrix njdm(dm.getSize());

  for(size_t i=0;i<dm.getSize();i++){
    Sequence_double data;
    data.dbl = -1;
    data.s.name = dm.getIdentifier(i);
    SequenceTree::Node *node = tree.getRoot()->addChild(data);
    njdm.setIdentifier(i,node);
    for(size_t j=i;j<dm.getSize();j++){
      njdm.setDistance(i,j,dm.getDistance(i,j));
    }
  }
*/
  //create NJ tree
  switch(m){
  case NJ:computeFloatNeighborJoiningTree(njflodm,defdata); break;
  case BIONJ: computeFloatBioNJTree(njflodm,defdata); break;
  case FNJ: computeFloatFNJTree(njflodm,defdata); break;
  default:
    PROG_ERROR("Unexpected method");
  }
}

