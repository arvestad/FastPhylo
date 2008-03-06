//--------------------------------------------------
//                                        
// File: TamuraNei.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: TamuraNei.cpp,v 1.3 2006/01/24 22:55:46 isaac Exp $                                 
//
//--------------------------------------------------

#include "dna_pairwise_sequence_likelihood.hpp"

#include <string>
#include <iostream>
#include <assert.h>

using namespace std;

//
// Computes the value of a parenthesis of the derivative of the
// likelihood.  Since we want to find the value that maximizes the
// likelihood we want to find the value x that find the zero of this
// expression.
//
// Here x represents B*t i.e. B the Markov-transition probability of a
// transversion multiplied with the time t.  Note that K = A/B where
// A is the Markov-transition probability of a transition. Thus if the
// transition/transversion ratio is 2 then K=4.


static float gA;
static float gC;
static float gG;
static float gT;
static float gU;
static float gY;

static float purine_ratio;//fix ratio
static float pyrimidine_ratio;//fix ratio
static float tU;//observed purine transitions
static float tY;//observed pyrimidine transitions
static float tV;//observed transversions
static float n; //the string length


static float E = exp(1.0);

// THIS IS ACTUALLY THE DERIVATIVE OF THE LOG LIKELIHOOD!!
static float
partof_derivative_likelihood(float b){
  float k1 = purine_ratio;
  float k2 = pyrimidine_ratio;
  float two_b = 2*b;
  float exp_2b = powf(E,two_b);
  
float VAL=
  (2*tV)/(-1 + exp_2b) +

  ///
  (((-2*(gY))/exp_2b + (2*(gY + (gU)*k1))/powf(E,two_b*(gY + (gU)*k1)))*tU)/
  (-powf(E,-two_b*(gY + (gU)*k1)) + gU + (gY)/exp_2b) + 
  ///
  (((-4*(gU)*(gY))/exp_2b -
    (2*gA*gG*((-2*(gY))/exp_2b + (2*(gY + (gU)*k1))/powf(E,two_b*(gY + (gU)*k1))))/(gU) -
    (2*gC*gT*((-2*(gU))/exp_2b + (2*(gU + (gY)*k2))/powf(E,two_b*(gU + (gY)*k2))))/(gY)
    )*(n - tV - tU - tY))/
  (1 - 2*(1 - powf(E,-two_b))*(gU)*(gY) -
   (2*gC*gT*(-powf(E,-two_b*(gU + (gY)*k2)) + gC + (gU)/exp_2b + gT))/(gY) - 
   (2*gA*gG*(-powf(E,-two_b*(gY + (gU)*k1)) + gU + (gY)/exp_2b))/(gU)) +

  ////
  (((-2*(gU))/exp_2b + (2*(gU + (gY)*k2))/powf(E,two_b*(gU + (gY)*k2)))*tY)/
  (-powf(E,-two_b*(gU + (gY)*k2)) + gC + (gU)/exp_2b + gT);

// cout << " VAL = " << VAL << endl;
 return VAL;
}
  
// static void
// TEST(){
//   cout << "------    TEST" << endl;
//   gA=0.25;
//   gC=0.25;
//   gG=0.25;
//   gT=0.25;
//   gU=0.5;
//   gY=0.5;
//   purine_ratio=4;
//   pyrimidine_ratio=4;
//   tU=2;
//   tY=2;
//   tV=2;
//   n=100;
//   //value computed with mathematica
//   cout << " PART OF " << partof_derivative_likelihood(0.9) << "      diff  " << fabs((-34.9386) - partof_derivative_likelihood(0.9))<<endl;
//   assert(fabs((-34.9386) - partof_derivative_likelihood(0.9))<=0.00005);
//   cout << "---" << endl;

//   gA=0.25;
//   gC=0.25;
//   gG=0.3;
//   gT=0.1;
//   gU=gA+gG;
//   gY=gC+gT;
//   purine_ratio=4;
//   pyrimidine_ratio=5;
//   tU=2;
//   tY=2;
//   tV=2;
//   n=100;

//   //value computed with mathematica
//   cout << " PART OF " << partof_derivative_likelihood(0.9) << "      diff  " << fabs((-19.4249) - partof_derivative_likelihood(0.9))<<endl;
//   assert(fabs((-19.4249) - partof_derivative_likelihood(0.9))<=0.00005);
//   cout << "------    ENDTEST" << endl;
// }


// static int maxiter=0;
// static int numabove =0;
// static int total=0;


// //Binary search for the zero of the derivative of the likelihood
// //which is the same as finding the zero of partof_derivative_likelihood.
// //This function is positive from [0->searchedval] and negative [searchedval->...].

// static float
// _binary_search(){
//   //initial guess
//   // A_u - markov transition probability of purines is approximated by the observed number of purine trans
//   // A_y - .... pyrimidine trans
//   // B  - markov transition prob of transversion... since there are two transversions B =~ tV/2
//   // note that A_u/B =
//   // thus to get a guess for B*t we use the following (A_u/k_u + A_y/k_y + B)/3
//   float bt_prob = ( (tU/purine_ratio) + (tY/pyrimidine_ratio) + tV)/(6.0*n);
//   //DALIG GISNING EFTERSOM DETTA INTE SAMMA BERAKNINGAR SOM I K2P.
//   float interval_length = bt_prob;
  
//   int i = 0;

//   if (  partof_derivative_likelihood(bt_prob) > 0 ){
//     do {
//       i++;
//       bt_prob += interval_length;
//       //cout << "initial " << bt_prob << "   ilen " << interval_length << endl;
//       //PENDING LOOP COUNTER!
//     } while ( partof_derivative_likelihood(bt_prob) > 0 );
//   }

//   //Now we know that the search value is between [bt_prob - interval_length, bt_prob]
//   // BEGIN BINARY SEARCH
//   while ( interval_length > 0.00001 ){
//     //PENDING LOOP COUNTER!
//     interval_length = interval_length*0.5;
//     float mid = bt_prob - interval_length;
//     if ( partof_derivative_likelihood(mid) < 0 )
//       bt_prob = mid;
//     i++;
//     //cout << "ITER " <<(i++) << "  "<< bt_prob << "   ilen " << interval_length << "    dist: " <<
//     //  4.0*(gA*gG*purine_ratio + gT*gC*pyrimidine_ratio + gU*gY)*bt_prob << endl;
//   }

//   if ( i > maxiter ){
//     maxiter = i;
//     cout << "max " << i << endl;
//   }
//   total++;
//   if ( i > 15 ){
//     numabove++;
//     cout << "iter " << i << "  " << ((float)numabove)/total << endl;
//   }
//   interval_length = interval_length*0.5;
//   bt_prob = bt_prob - interval_length;


//   //cout << "binary search     val: " << bt_prob << "     funcval  " << partof_derivative_likelihood(bt_prob) << endl;
  
//   return bt_prob;
// }



static float
_secant_search(){
  float x0 = ( (tU/purine_ratio) + (tY/pyrimidine_ratio) + tV)/(6.0*n);
  float fx0 = partof_derivative_likelihood(x0);
  float x1;
  if ( fx0 < 0 )
    x1 = x0/2;
  else
    x1 = x0*3/2;
  
  //cout << "fx0 " << fx0 << "   x0 " << x0 << endl; 
  int i =0;
  while (i < 20 && fabs(x1-x0) > 0.00001){
    float tmp = x1;
    float tmpf = partof_derivative_likelihood(x1);

    x1 = x1 - (x1-x0)/(tmpf-fx0)*tmpf;
    fx0= tmpf;
    x0 = tmp;
    //cout << "fx0 " << fx0 << "   x0 " << x0 << endl;
    i++;
  }
  
  
//   if ( i > maxiter ){
//     maxiter = i;
//     cout << "max " << i << endl;
//   }
//   total++;
//   if ( i > 15 ){
//     numabove++;
//     cout << "iter " << i << "  " << ((float)numabove)/total << endl;
//   }

  //cout << "secant search     val: " << x1 << "     funcval  " << partof_derivative_likelihood(x1) << endl; 
  return x1;
}



ML_string_distance
compute_Tamura_Nei_fixratio(int strlen, TN_string_distance sd,
                            int numAs, int numCs,
                            int numGs, int numTs,
                            float purine_ts_tv_ratio,
                            float pyrimidine_ts_tv_ratio){
  purine_ratio = 2.0* purine_ts_tv_ratio; //K is A/B. while fixRatio is A/2B
  pyrimidine_ratio = 2.0* pyrimidine_ts_tv_ratio; //K is A/B. while fixRatio is A/2B
  tU = sd.purine_transitions;
  tY = sd.pyrimidine_transitions;
  tV = sd.transversions;
  n = (1.0*strlen) - sd.deletedPositions;
  
  float norm = numAs + numCs + numGs + numTs;
  
  gA = ((float)numAs)/norm;
  gA = ( gA < 0.000001 ? 0.000001 : gA);
  gC = 0.000001+((float)numCs)/norm;
  gC = ( gC < 0.000001 ? 0.000001 : gC);
  gG = 0.000001+((float)numGs)/norm;
  gG = ( gG < 0.000001 ? 0.000001 : gG);
  gT = 0.000001+((float)numTs)/norm;
  gT = ( gT < 0.000001 ? 0.000001 : gT);
  
  gU = gA+gG;
  gY = gC+gT;

  //cout << "gA " << gA << endl;
  //cout << "gC " << gC << endl;
  //cout << "gG " << gG << endl;
  //cout << "gT " << gT << endl;
  
  float bt_prob;
  //bt_prob = _binary_search();
  bt_prob = _secant_search();
  
  

  //the distance
  ML_string_distance tp;
  tp.distance = 4.0*(gA*gG*purine_ratio + gT*gC*pyrimidine_ratio + gU*gY)*bt_prob;

  // FIX THE CHANGE PROBABILITIES

  float tv_prob = 2.0*gU*gY*(1-exp(-2*bt_prob));
  float ts_pyrimidine_prob = 2.0 * gT*gC/gY*(gY- exp(-2*(gY*pyrimidine_ratio+ gU)*bt_prob) + gU*exp(-2*bt_prob));
  float ts_purine_prob = 2.0 * gA*gG/gU*(gU- exp(-2*(gU*purine_ratio+ gY)*bt_prob) + gY*exp(-2*bt_prob));

  assert ( tv_prob > 0.00001 && tv_prob < 0.99999);
  assert (ts_pyrimidine_prob > 0.00001 &&  ts_pyrimidine_prob < 0.99999);
  assert ( ts_purine_prob > 0.00001 && ts_purine_prob < 0.99999);
  
  tp.A_A = 1.0-ts_purine_prob-tv_prob; tp.A_C = tv_prob*0.5;  tp.A_G = ts_purine_prob; tp.A_T = tv_prob*0.5;
  tp.C_C = 1.0-ts_pyrimidine_prob-tv_prob;  tp.C_G = tv_prob*0.5; tp.C_T = ts_pyrimidine_prob;
  tp.G_G = 1.0-ts_purine_prob-tv_prob;  tp.G_T = tv_prob*0.5;
  tp.T_T = 1.0-ts_pyrimidine_prob-tv_prob;

  //  TEST();
 
  return tp;
  
}
