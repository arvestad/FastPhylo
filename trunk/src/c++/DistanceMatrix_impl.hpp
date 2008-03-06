//--------------------------------------------------
//                                        
// File: DistanceMatrix_impl.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: DistanceMatrix_impl.hpp,v 1.10 2006/12/19 12:36:23 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef DISTANCEMATRIX_IMPL_HPP
#define DISTANCEMATRIX_IMPL_HPP

#include <string>
#include "file_utils.hpp"

DM_TEMPLATE void
DISTANCEMATRIX::assureSize(){
  identifiers.resize(size);
  D.resize(size);
  for (size_t i = 0 ; i < size ; i++ )
    D[i].resize(size);
}

DM_TEMPLATE
DISTANCEMATRIX::DistanceMatrix(size_t size) {
  this->size = size;
  assureSize();
}


DM_TEMPLATE
DISTANCEMATRIX::DistanceMatrix(const DISTANCEMATRIX &dm){
  (*this) = dm;
}

DM_TEMPLATE DISTANCEMATRIX&
DISTANCEMATRIX::operator=(const DISTANCEMATRIX&dm){
  size = dm.size;

  identifiers.resize(size);
  D.resize(size);

  //copy the upper triangle
  for ( size_t i = 0 ; i < size ; i++ ){
    identifiers[i] = dm.identifiers[i];
    DISTVEC &thisvec = D[i];
    thisvec.resize(size);
    const DISTVEC &thatvec = dm.D[i];
    for ( size_t j = i ; j < size ; j++ )
      thisvec[j] = thatvec[j];
  }

  return *this;
}


DM_TEMPLATE
DISTANCEMATRIX::DistanceMatrix(std::istream &in){
  size = -1;
  objInitFromStream(in);
}
  

DM_TEMPLATE void
DISTANCEMATRIX::setDefaultValues(DistanceType &defval, Identifier &defid){

  for ( size_t i = 0 ; i < size ; i++ )
    identifiers[i] = defid;

  for ( size_t i = 0 ; i < size ; i++ ){
    for ( size_t j = 0 ; j <= i ; j++ )
      D[j][i] = defval;
  }
}


//CREATE FROM STREAM
// here is an example of a phylip file
// Should be given on full form
// the diagonal should be all 0s
//
//     7
// Bovine      0.0000  1.6866  1.7198  1.6606  1.5243  1.6043  1.5905
// Mouse       1.6866  0.0000  1.5232  1.4841  1.4465  1.4389  1.4629
// Gibbon      1.7198  1.5232  0.0000  0.7115  0.5958  0.6179  0.5583
// Orang       1.6606  1.4841  0.7115  0.0000  0.4631  0.5061  0.4710
// Gorilla     1.5243  1.4465  0.5958  0.4631  0.0000  0.3484  0.3083
// Chimp       1.6043  1.4389  0.6179  0.5061  0.3484  0.0000  0.2692
// Human       1.5905  1.4629  0.5583  0.4710  0.3083  0.2692  0.0000


DM_TEMPLATE std::istream&
DISTANCEMATRIX::objInitFromStream(std::istream &in){

  size_t newSize;
  in >> newSize;

  if ( newSize != size ){
    size = newSize;
    assureSize();
  }
  
  // read each line of the matrix and set the distances
  for ( size_t i = 0 ; i < size ; i++ ){
    idInit(in,identifiers[i]);

    DistanceType dist;
    
    //skip the lower left part
    size_t j = 0;
    while ( j < i ){
      j++;
      distInit(in,dist);
    }

    distInit(in,D[i][i]);

    j++;
    for ( ; j < size ; j++ ){
      distInit(in,D[i][j]);
    }
  }

  return in;
}


DM_TEMPLATE std::ostream&
DISTANCEMATRIX::printOn(std::ostream &out) const{
  out << std::setw(5) << std::right << size << std::endl;

  out.precision(6);
  for ( size_t i = 0 ; i < size ; i++ ){
    out  << std::setw(10) << std::left;
    idPrintOn(out,identifiers[i]);
    out << "  ";
    size_t j = 0;
    for (  ; j < i ; j++ ){
      out  << std::setw(10) << std::right;
      distPrintOn(out,D[j][i]);
      out << " ";
    }
    for ( ; j < size ; j++ ){
      out  << std::setw(10) << std::right;
      distPrintOn(out,D[i][j]);
      out << " ";
    }
    out << std::endl;
  }

  return out;
}


//SWAP AND REMOVE LAST ROW
DM_TEMPLATE void
DISTANCEMATRIX::swapRowToLast(size_t row){

  assert( row < size );
  size_t lastRow = size-1;
  if ( row == lastRow )
    return;

  //switch the distances
  for ( size_t i = 0 ; i < row ; i++ ){
    DistanceType tmp = D[i][row];
    D[i][row] = D[i][lastRow];
    D[i][lastRow] = tmp;
  }

  DistanceType tmp = D[row][row];
  D[row][row] = D[lastRow][lastRow];
  D[lastRow][lastRow] = tmp;

  for ( size_t i = row+1 ; i < lastRow ; i++ ){
    DistanceType tmp = D[row][i];
    D[row][i] = D[i][lastRow];
    D[i][lastRow] = tmp;
  }
  
  // Switch the Identifiers
  Identifier tmpI = identifiers[row];
  identifiers[row] = identifiers[lastRow];
  identifiers[lastRow] = tmpI;  
}

DM_TEMPLATE void
DISTANCEMATRIX::removeLastRow(){
  size--;
}




#endif // DISTANCEMATRIX_IMPL_HPP
