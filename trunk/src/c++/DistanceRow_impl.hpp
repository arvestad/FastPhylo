//--------------------------------------------------
//                                        
// File: DistanceRow_impl.hpp                             
//                             
// Author: Richard Schobesberger         
// e-mail: rschob@sbs.su.se   
//                             
//                               
//
//--------------------------------------------------
#ifndef DISTANCEROW_IMPL_HPP
#define DISTANCEROW_IMPL_HPP

#include <string>
#include "file_utils.hpp"

DM_TEMPLATE void
DISTANCEROW::assureSize(){
	//Immer eine Zeile
  D.resize(columns);
}

DM_TEMPLATE
DISTANCEROW::DistanceRow(size_t columns) {
  this->columns = columns;
  assureSize();
}


DM_TEMPLATE
DISTANCEROW::DistanceRow(const DISTANCEROW &dm){
  (*this) = dm;
}


//Doesn't function here :)
DM_TEMPLATE
DISTANCEROW::DistanceRow(std::istream &in){
  columns = -1;
  objInitFromStream(in);
}
  
//Doesn't function here :)
DM_TEMPLATE void
DISTANCEROW::setDefaultValues(DistanceType &defval, Identifier &defid){

	identifier = defid;

	for ( size_t i = 0 ; i < columns ; i++ ){
		D[i] = defval;
	}
}


#endif // DISTANCEROW_IMPL_HPP
