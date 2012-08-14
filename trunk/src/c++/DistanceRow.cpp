#include <float.h>
#include <math.h>
#include <iostream>
#include "log_utils.hpp"
#include "DistanceRow.hpp"
#include "config.h"
#include <string>

using std::string;
using std::cout;
using std::endl;

bool
applyFixFactorRow(StrFloRow &dm, float fixFactor){
  float biggest = 0;
  int size =dm.getColumns();

	for ( int j = 0 ; j < size ; j++ ){
	  float d = dm.getDistance(j);
	  if ( isfinite(d) && d>biggest ){
	    biggest = d;
	  }
  	}

  bool changed = false;
  biggest = biggest*fixFactor;
	for ( int j = 0 ; j < size ; j++ ){
	  float d = dm.getDistance(j);
	  if ( !isfinite(d) || d<0 ){
	    dm.setDistance(j,biggest);
		changed  = true;
	  }
  	}

  return changed;
}  
