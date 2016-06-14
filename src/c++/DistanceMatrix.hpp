//--------------------------------------------------
//                                        
// File: DistanceMatrix.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: DistanceMatrix.hpp,v 1.12 2006/12/08 11:09:12 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef DISTANCEMATRIX_HPP
#define DISTANCEMATRIX_HPP

#include "Object.hpp"
#include <vector>
#include <string>
#include "InitAndPrintOn_utils.hpp"

//
// A symetric distance matrix. There is one template parameter
// describing the data in the elements of the matrix and one template
// parameter describing the data used as identifier for each
// row/column. In addition, these types have two be provided with
// function objects for initiating and printing the data types to a
// stream.
// Ex. DistanceMatrix<std::string, double, 
//                   Data_init<std::string>, Data_printOn<std::string>, 
//                   Data_init<double>, Data_printOn<double> >
//
//
#define DISTANCEMATRIX DistanceMatrix<Identifier,DistanceType,IdentInit,IdentPrintOn,DistInit,DistPrintOn>
#define DM_TEMPLATE template<class Identifier, class DistanceType, class IdentInit, class IdentPrintOn, class DistInit, class DistPrintOn>

#define DISTVEC std::vector<DistanceType>

// A SYMMETRIC DISTANCE MATRIX

DM_TEMPLATE
class DistanceMatrix : public Object{
public:

  //CONSTRUCTORS
  DistanceMatrix(size_t size=0);
  DistanceMatrix(const DistanceMatrix &dm);
  DistanceMatrix& operator=(const DistanceMatrix &dm);
  DistanceMatrix(std::istream &in);

  //DIMENSIONS
  size_t getSize() const{ return size;}
  void resize(size_t newsize) { if ( size != newsize) {size = newsize; assureSize();} }
  
  //fills it from the stream. The stream should be on a
  //regular phylip format.
  virtual std::istream& objInitFromStream(std::istream &is);

  
  //set all the matrix elements and identifiers to the supplied values.
  void setDefaultValues(DistanceType &defval, Identifier &defid);
  
  //---------------------------------
  //GET AND SET DISTANCE
  DistanceType getDistance(int i, int j) const{
    //only the upper right triangle 
    if ( i <= j )
      return D[i][j];
    else
      return D[j][i];
  };
  
  void setDistance(int i, int j, DistanceType d) {
    if ( i <= j )
      D[i][j] = d;
    else
      D[j][i] = d;
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
  
  //----------------------
  std::ostream& printOn(std::ostream &out) const;

  //------------------- PRIVATE ------------------------
private:
  DistInit distInit;
  DistPrintOn distPrintOn;
  IdentInit idInit;
  IdentPrintOn idPrintOn;
  
  size_t size;
  std::vector<Identifier> identifiers;
  std::vector<std::vector<DistanceType> > D;

  //makes sure that the vectors are of the right size.
  void assureSize();

};


//--------------------------
// Include the implementation
#include "DistanceMatrix_impl.hpp"


//-------------------------
// A regular distance matrix with strings as identifier and doubles as distance type.
typedef DistanceMatrix<std::string, double, Data_init<std::string>, Data_printOn<std::string>, Data_init<double>, Data_printOn<double> > StrDblMatrix;

//Distance for over saturated data ( !isfinite(d) || d<0 ) is set to <int>*maxRealDistance.\n"
//returns true if some element was changed. 
bool
applyFixFactor(StrDblMatrix &dm, double fixFactor);


#endif // DISTANCEMATRIX_HPP






