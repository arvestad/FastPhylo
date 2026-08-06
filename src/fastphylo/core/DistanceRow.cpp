#include <cfloat>
#include <cmath>
#include <iostream>
#include "fastphylo/core/log_utils.hpp"
#include "fastphylo/core/DistanceRow.hpp"
#include "config.h"
#include <string>

using std::string;
using std::cout;
using std::endl;

bool
applyFixFactorRow(StrFloRow &dm, float fixFactor){
  float biggest = 0;
  size_t size =dm.getColumns();

	for ( size_t j = 0 ; j < size ; j++ ){
	  float d = dm.getDistance(j);
	  if ( std::isfinite(d) && d>biggest ){
	    biggest = d;
	  }
  	}

  bool changed = false;
  biggest = biggest*fixFactor;
	for ( size_t j = 0 ; j < size ; j++ ){
	  float d = dm.getDistance(j);
	  if ( !std::isfinite(d) || d<0 ){
	    dm.setDistance(j,biggest);
		changed  = true;
	  }
  	}

  return changed;
}  
