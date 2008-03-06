
#include "std_c_utils.h"

#include <stdlib.h>

int
randomInt( int a, int b){
  return a+(int) (b*1.0*rand()/(RAND_MAX+1.0));
}


float
randomFloat(float a, float b) {
  return a + (((float)rand())/RAND_MAX) * ( b - a );
}


