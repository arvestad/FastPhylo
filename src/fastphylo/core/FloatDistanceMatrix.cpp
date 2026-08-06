/*
 * FloatDistanceMatrix.cpp
 *
 *  Created on: Dec 7, 2011
 *      Author: Mehmood Alam Khan
 *      Email : malagori@kth.se
 */

#include "fastphylo/core/FloatDistanceMatrix.hpp"

#include <cfloat>
#include <cmath>
#include <iostream>
#include "fastphylo/core/log_utils.hpp"
#include "config.h"
#include <string>

using std::string;
using std::cout;
using std::endl;

bool
applyFixFactor(StrFloMatrix &dm, float fixFactor){
  float biggest = 0;
  int size =dm.getSize();

  for ( int i = 0 ; i < size ; i++ ){
    for ( int j = 0 ; j < size ; j++ ){
      float d = dm.getDistance(i,j);
      if ( std::isfinite(d) && d>biggest ){
        biggest = d;
      }
    }
  }

  bool changed = false;
  biggest = biggest*fixFactor;
  for ( int i = 0 ; i < size ; i++ ){
    for ( int j = 0 ; j < size ; j++ ){
      float d = dm.getDistance(i,j);
      if ( !std::isfinite(d) || d<0 ){
        dm.setDistance(i,j,biggest);
	changed  = true;
      }
    }
  }

  return changed;
}

