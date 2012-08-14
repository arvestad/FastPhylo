/*
 * FloatDistanceMatrix.hpp
 *
 *  Created on: Dec 7, 2011
 *      Author: Mehmood Alam Khan
 *      Email : malagori@kth.se
 */
#ifndef FLOATDISTANCEMATRIX_IMPL_HPP
#define FLOATDISTANCEMATRIX_IMPL_HPP

#include <string>
#include "file_utils.hpp"

DM_TEMPLATE void
FLOATDISTANCEMATRIX::assureSize(){
  identifiers.resize(rows);
  D.resize(rows);

  for (size_t i = 0 ; i < rows ; i++ ) {

  	if(D[i] == 0) {
  	  D[i] = new std::vector<DistanceType>;
  	}

  	(*D[i]).resize(columns - i);

  }
}

DM_TEMPLATE
FLOATDISTANCEMATRIX::FloatDistanceMatrix(size_t columns) {
  this->columns = columns;
  this->rows = columns;
  assureSize();
}

DM_TEMPLATE
FLOATDISTANCEMATRIX::FloatDistanceMatrix(size_t rows, size_t columns) {
  this->columns = columns;
  this->rows = rows;
  assureSize();
}

DM_TEMPLATE
FLOATDISTANCEMATRIX::FloatDistanceMatrix(const FLOATDISTANCEMATRIX &dm){
  (*this) = dm;
}

DM_TEMPLATE FLOATDISTANCEMATRIX&
FLOATDISTANCEMATRIX::operator=(const FLOATDISTANCEMATRIX &dm){

  for(int i = 0; i < rows; ++i) {
  	 if(D[i] != 0) {
  	 	delete D[i];
  	 }
  }

  rows = dm.rows;
  columns = dm.columns;

  identifiers.resize(rows);
  D.resize(rows);

  //copy the upper triangle
  for ( size_t i = 0 ; i < rows ; i++ ){
    identifiers[i] = dm.identifiers[i];

    DISTVEC &thisvec = (*D[i]);

    thisvec.resize(i);
    const DISTVEC &thatvec = (*dm.D[i]);

    for ( size_t j = 0 ; j < rows - i ; j++ ) {
      thisvec[j] = thatvec[j];
    }
  }

  return *this;
}


DM_TEMPLATE
FLOATDISTANCEMATRIX::FloatDistanceMatrix(std::istream &in){
  rows = -1;
  columns = -1;
  objInitFromStream(in);
}


DM_TEMPLATE void
FLOATDISTANCEMATRIX::setDefaultValues(DistanceType &defval, Identifier &defid){

  for ( size_t i = 0 ; i < rows ; i++ ) {
    identifiers[i] = defid;
  }

  for ( size_t i = 0 ; i < rows ; i++ ){
    for ( size_t j = 0 ; j < columns - i; j++ ) {
      (*D[i])[j] = defval;
    }
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
FLOATDISTANCEMATRIX::objInitFromStream(std::istream &in){

  size_t newSize;
  in >> newSize;

  if ( newSize != rows || newSize != columns ){
    rows = newSize;
    columns = newSize;
    assureSize();
  }

  // read each line of the matrix and set the distances
  for ( size_t i = 0 ; i < rows ; i++ ){
    idInit(in,identifiers[i]);

    DistanceType dist;

    //skip the lower left part
    size_t j = 0;
    while ( j < i ){
      j++;
      distInit(in, dist);
    }

    distInit(in, ((*D[i])[i - i]));

    j++;
    for (; j < columns ; j++ ){
      distInit(in, ((*D[i])[j - i]));
    }
  }
  return in;
}

DM_TEMPLATE std::ostream&
FLOATDISTANCEMATRIX::printOn(std::ostream &out) const{

  out << std::setw(5) << std::right << rows << std::endl;

  out.precision(6);
  for ( size_t i = 0 ; i < rows ; i++ ){
    out  << std::setw(10) << std::left;
    idPrintOn(out,identifiers[i]);
    out << "  ";
    size_t j = 0;
    for (  ; j < i ; j++ ){
      out  << std::setw(10) << std::right;
      distPrintOn(out,(*D[j])[i - j]);
      out << " ";
    }
    for ( ; j < columns ; j++ ){
      out  << std::setw(10) << std::right;
      distPrintOn(out,(*D[i])[j - i]);
      out << " ";
    }
    out << std::endl;
  }
  return out;
}


//SWAP AND REMOVE LAST ROW
DM_TEMPLATE void
FLOATDISTANCEMATRIX::swapRowToLast(size_t row){

  assert( row < rows );
  size_t lastRow = rows - 1;
  if ( row == lastRow ) {
    return;
  }

  //switch the distances

  if(lastRow < row) {
  	(*D[row]).resize(columns - lastRow);
  } else {
  	(*D[lastRow]).resize(columns - row);
  }

  //change the values in every line before the line which shall be moved
  for ( size_t i = 0 ; i < row ; i++ ){
    DistanceType tmp = (*D[i])[row - i];
    (*D[i])[row - i] = (*D[i])[lastRow - i];
    (*D[i])[lastRow - i] = tmp;
  }

  //change the zeros
  DistanceType tmp = (*D[row])[row - row];
  (*D[row])[row - row] = (*D[lastRow])[lastRow - lastRow];
  (*D[lastRow])[lastRow - lastRow] = tmp;

  //change the values in the columns
  for ( size_t i = row+1 ; i < lastRow ; i++ ){
    DistanceType tmp = (*D[row])[i - row];
    (*D[row])[i - row] = (*D[i])[lastRow - i];
    (*D[i])[lastRow - i] = tmp;
  }

  // Switch the Identifiers
  Identifier tmpI = identifiers[row];
  identifiers[row] = identifiers[lastRow];
  identifiers[lastRow] = tmpI;
}

DM_TEMPLATE void
FLOATDISTANCEMATRIX::removeLastRow(){
  if(D[rows - 1] != 0) {
  	 delete D[rows - 1];
  }
  --rows;
  assureSize();
}


DM_TEMPLATE void
FLOATDISTANCEMATRIX::addRow(){
	if(rows < columns){
		++rows;
	}
	assureSize();
}


DM_TEMPLATE void
FLOATDISTANCEMATRIX::makeQuadratic(){}
#endif // FLOATDISTANCEMATRIX_IMPL_HPP
