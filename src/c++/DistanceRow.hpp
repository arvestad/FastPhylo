//--------------------------------------------------
//                                        
// File: DistanceRow.hpp                             
//                             
// Author: Richard Schobesberger         
// e-mail: rschob@sbc.su.se
//                             
//                            
//
//--------------------------------------------------
#ifndef DISTANCEROW_HPP
#define DISTANCEROW_HPP

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
#define DISTANCEROW DistanceRow<Identifier,DistanceType,IdentInit,IdentPrintOn,DistInit,DistPrintOn>
#define DM_TEMPLATE template<class Identifier, class DistanceType, class IdentInit, class IdentPrintOn, class DistInit, class DistPrintOn>

#define DISTVEC std::vector<DistanceType>

// A SYMMETRIC DISTANCE MATRIX

DM_TEMPLATE
class DistanceRow : public Object{
public:

  //CONSTRUCTORS
  DistanceRow(size_t columns=0);
  DistanceRow(const DistanceRow &dm);
  DistanceRow(std::istream &in);

  //DIMENSIONS
  size_t getColumns() const{ return columns;}
  void resize(size_t newsize) { if ( columns != newsize) {columns = newsize; assureSize();} }
  
  //set all the matrix elements and identifiers to the supplied values.
  void setDefaultValues(DistanceType &defval, Identifier &defid);
  
  //---------------------------------
  //GET AND SET DISTANCE
  DistanceType getDistance(int j) const{
      return D[j];
  };
  
  void setDistance(int j, DistanceType d) {
      D[j] = d;
  };

  //---------------------
  void setIdentifier(Identifier id){
    identifier = id;
  }
  Identifier& getIdentifier(){
    return identifier;
  }
 const Identifier& getIdentifier() const{
    return identifier;
 }

  //------------------- PRIVATE ------------------------
private:
  DistInit distInit;
  DistPrintOn distPrintOn;
  IdentInit idInit;
  IdentPrintOn idPrintOn;
  
  size_t columns;
  Identifier identifier;
  std::vector<DistanceType> D;

  //makes sure that the vectors are of the right size.
  void assureSize();

};


//--------------------------
// Include the implementation
#include "DistanceRow_impl.hpp"


//-------------------------
// A regular distance matrix with strings as identifier and doubles as distance type.
typedef DistanceRow<std::string, double, Data_init<std::string>, Data_printOn<std::string>, Data_init<double>, Data_printOn<double> > StrDblRow;
typedef DistanceRow<std::string, float, Data_init<std::string>, Data_printOn<std::string>, Data_init<float>, Data_printOn<float> > StrFloRow;
//Distance for over saturated data ( !isfinite(d) || d<0 ) is set to <int>*maxRealDistance.\n"
//returns true if some element was changed. 
bool
applyFixFactorRow(StrFloRow &dm, float fixFactor);



#endif // DISTANCEROW_HPP






