#include <float.h>
#include <math.h>
#include "NeighborJoining.hpp"
#include "Sequence.hpp"

typedef DistanceMatrix<SequenceTree::Node *,double,Data_init<SequenceTree::Node *>,Data_printOn<SequenceTree::Node *>,Data_init<double>,Data_printOn<double> > NJMatrix;


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


