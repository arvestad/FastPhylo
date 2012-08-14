/*
 * FloatDistanceMatrix.hpp
 *
 *  Created on: Dec 7, 2011
 *      Author: Mehmood Alam Khan
 *      Email : malagori@kth.se
 */

#ifndef FLOATDISTANCEMATRIX_HPP_
#define FLOATDISTANCEMATRIX_HPP_

#include "Object.hpp"
#include <vector>
#include <string>
#include "InitAndPrintOn_utils.hpp"

//
// An A-symetric distance matrix. There is one template parameter
// describing the data in the elements of the matrix and one template
// parameter describing the data used as identifier for each
// row/column. In addition, these types have two be provided with
// function objects for initiating and printing the data types to a
// stream.
// Ex. FloatDistanceMatrix<std::string, double,
//                   Data_init<std::string>, Data_printOn<std::string>,
//                   Data_init<double>, Data_printOn<double> >
//
//
#define FLOATDISTANCEMATRIX FloatDistanceMatrix<Identifier,DistanceType,IdentInit,IdentPrintOn,DistInit,DistPrintOn>
#define DM_TEMPLATE template<class Identifier, class DistanceType, class IdentInit, class IdentPrintOn, class DistInit, class DistPrintOn>

#define DISTVEC std::vector<DistanceType>

// A SYMMETRIC DISTANCE MATRIX

DM_TEMPLATE
class FloatDistanceMatrix : public Object{
public:

  //CONSTRUCTORS
  FloatDistanceMatrix(size_t columns=0);
  FloatDistanceMatrix(size_t rows, size_t columns);
  FloatDistanceMatrix(const FloatDistanceMatrix &dm);
  FloatDistanceMatrix& operator=(const FloatDistanceMatrix &dm);
  FloatDistanceMatrix(std::istream &in);

  //DESTRUCTOR
  ~FloatDistanceMatrix() {
  	 for(int i = 0; i < rows; ++i) {
	   if(D[i] != 0) {
	     delete D[i];
	   }
  	 }
  }

  //DIMENSIONS
  size_t getSize() const{ return rows;}
  size_t getRows() const{ return rows;}
  size_t getColumns() const{return columns;}

  void resize(size_t newRows, size_t newColumns) {
  	if (rows != newRows || columns != newColumns) {
  		rows = newRows;
  		columns = newColumns;
  		assureSize();
  	}
  }

  void resize(size_t newSize) {
  	resize(newSize, newSize);
  }

  //fills it from the stream. The stream should be on a
  //regular phylip format.
  virtual std::istream& objInitFromStream(std::istream &is);

  //fills it from the stream which is binary
  //for format details check BinaryOutputStream.hpp/cpp in programs/fastdist
  //virtual void objInitFromBinaryStream(std::istream *in);


  //set all the matrix elements and identifiers to the supplied values.
  void setDefaultValues(DistanceType &defval, Identifier &defid);

  //---------------------------------
  //GET AND SET DISTANCE
  DistanceType getDistance(int i, int j) const{
    //only the upper right triangle
    if ( i <= j ) {
      return (*D[i])[j - i];
    } else {
      return (*D[j])[i - j];
    }
  };

  void setDistance(int i, int j, DistanceType d) {
    if ( i <= j ) {
    	(*D[i])[j - i] = d;
    } else {
      (*D[j])[i - j] = d;
    }
  };

  //---------------------
  void setIdentifiers(std::vector<Identifier> ids){
    for(size_t i=0;i<ids.size();i++)
      identifiers[i]=ids[i];
  }
  void setIdentifier(int i, Identifier id){
    identifiers[i] = id;
  }
  Identifier& getIdentifier(int i){
    return identifiers[i];
  }
 const Identifier& getIdentifier(int i) const{
    return identifiers[i];
 }

  //SWAP AND REMOVE LAST ROW
  void swapRowToLast(size_t row);
  void removeLastRow();
  void addRow();
  void makeQuadratic();

  //----------------------
  std::ostream& printOn(std::ostream &out) const;

  //------------------- PRIVATE ------------------------
private:
  DistInit distInit;
  DistPrintOn distPrintOn;
  IdentInit idInit;
  IdentPrintOn idPrintOn;

  size_t columns;
  size_t rows;
  std::vector<Identifier> identifiers;
  std::vector<std::vector<DistanceType>*> D;

  //makes sure that the vectors are of the right size.
  void assureSize();

};


//--------------------------
// Include the implementation
#include "FloatDistanceMatrix_impl.hpp"


//-------------------------
// A regular distance matrix with strings as identifier and doubles as distance type.
//typedef FloatDistanceMatrix<std::string, double, Data_init<std::string>, Data_printOn<std::string>, Data_init<double>, Data_printOn<double> > StrDblMatrix;
typedef FloatDistanceMatrix<std::string, float, Data_init<std::string>, Data_printOn<std::string>, Data_init<float>, Data_printOn<float> > StrFloMatrix;

//Distance for over saturated data ( !isfinite(d) || d<0 ) is set to <int>*maxRealDistance.\n"
//returns true if some element was changed.
bool
applyFixFactor(StrFloMatrix &dm, float fixFactor);


#endif /* FLOATDISTANCEMATRIX_HPP_ */
