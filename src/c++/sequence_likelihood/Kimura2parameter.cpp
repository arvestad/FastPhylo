//--------------------------------------------------
//                                        
// File: Kimura2Parameter.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: Kimura2parameter.cpp,v 1.4 2006/12/10 19:58:53 isaac Exp $
//
//--------------------------------------------------

#include "dna_pairwise_sequence_likelihood.hpp"
#include <string>
#include <iostream>
#include <assert.h>

using namespace std;

//
// This file contains the function compute_K2P_fixratio() which
// computes the most likely estimate for the evolutionary distance
// under the Kimura 2 parameter model.
//




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
static float
partof_derivative_likelihood(float x,
                             float K,//fix ratio
                             float a,//observed transitions
                             float b,//observed transversions
                             float n){ //the string length


  float exp2x = exp(2.0*x);
  float exp4x = exp2x*exp2x;
  float exp2Kx =  exp(2.0*x*K);


  float VAL = //
    (b/(exp4x-1.0)) //
    + ( a*( exp2x*(1.0+K) - exp2Kx) / ( exp2Kx*(1.0+exp4x) - 2.0*exp2x) ) //
    + ( ( a+b-n)*( exp2x*(1.0+K)+exp2Kx ) / ( exp2Kx*(1.0+exp4x) + 2.0*exp2x) );

  //  cout << " VAL = " << VAL << endl;
  return VAL;
}


static float
_secant_search(
              float K,//fix ratio
              float a,//observed transitions
              float b,//observed transversions
              float n){
  float x0 = (a/K + b/2)/(n);
  float fx0 = partof_derivative_likelihood(x0,K,a,b,n);
  float x1;
  if ( fx0 < 0 )
    x1 = x0/2;
  else
    x1 = x0*3/2;
  
  //cout << "fx0 " << fx0 << "   x0 " << x0 << endl; 
  int i =0;
  while (i < 20 && fabs(x1-x0) > 0.00001){
    float tmp = x1;
    float tmpf = partof_derivative_likelihood(x1,K,a,b,n);

    x1 = x1 - (x1-x0)/(tmpf-fx0)*tmpf;
    assert ( x1 >= 0 );
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

  return x1;
}





ML_string_distance
compute_K2P_fixratio(int strlen, simple_string_distance sd, float fixRatio){
  
  float K = 2.0* fixRatio; //K is A/B. while fixRatio is A/2B
  float a = sd.transitions;
  float b = sd.transversions;

  float n = (1.0*strlen) - sd.deletedPositions;

  ML_string_distance tp;
  float bt_prob;
  float at_prob;
  
  //if the number of transitions or transversions is 0 then we can not compute
  //using the fix ration. Instead we just return the regular K2P distance.
  float a_div_b = a/b;
  if ( a_div_b > 10000 || a_div_b < 0.00001 ){
    bt_prob = b/n;//PENDING we should try using K here
    at_prob = a/n;
    tp.distance = -0.5 * log(1.0 - 2*bt_prob - at_prob) - 0.25 *log(1 - 2.0*at_prob);
  }
  else {
    //float bt_prob = _binary_search(K,a,b,n);
    bt_prob = _secant_search(K,a,b,n);//about 45% of the computation time is spent here
    //float bt_prob = _NEWTON_RAPHSON(K,a,b,n);
    //the distance
    tp.distance =   (K+2)*bt_prob;
    at_prob = bt_prob*K;
  }

  // FIX THE CHANGE PROBABILITIES
  assert(at_prob >= 0.0);
  assert(bt_prob >= 0.0);

  float tv_prob = 0.5*(1-exp(-4*bt_prob));
  float ts_prob = 0.25*(1-2*exp(-2*(at_prob+bt_prob)) + exp(-4*bt_prob));
  
  float id_prob = 1 - tv_prob - ts_prob;
  assert ( tv_prob >= 0.0 && tv_prob <= 1.0);
  assert ( ts_prob >= 0.0 && ts_prob <= 1.0);
  assert ( id_prob >= 0.0 && id_prob <= 1.0);
  
  tp.A_A = id_prob; tp.A_C = tv_prob*0.5;  tp.A_G = ts_prob; tp.A_T = tv_prob*0.5;
  tp.C_C = id_prob;  tp.C_G = tv_prob*0.5; tp.C_T = ts_prob;
  tp.G_G = id_prob;  tp.G_T = tv_prob*0.5;
  tp.T_T = id_prob;
  return tp;
  
}




// static float
// SLOW_partof_derivative_likelihood(float x,
//                                   float K,//fix ratio
//                                   float a,//observed transitions
//                                   float b,//observed transversions
//                                   float n){ //the string length


//   float E = exp(1);

//   return b/(-1 + powf(E,4*x)) + (a*(-powf(E,2*(1 + K)*x) + powf(E,4*x)*(1 + K)))/
//     (-2*powf(E,4*x) + powf(E,2*(1 + K)*x) + powf(E,2*(3 + K)*x)) + 
//    ((powf(E,2*(1 + K)*x) + powf(E,4*x)*(1 + K))*(a + b - n))/
//     (2*powf(E,4*x) + powf(E,2*(1 + K)*x) + powf(E,2*(3 + K)*x));
// }

// static float
// NEWTON_RAPHSON_RATIO(float x,
//                      float K,//fix ratio
//                      float a,//observed transitions
//                      float b,//observed transversions
//                      float n){ //the string length

//   float exp2x = exp(2.0*x);
//   float exp4x = exp2x*exp2x;
//   float exp2Kx =  exp(2.0*x*K);

  
//   float VAL = //
//     (b/(exp4x-1.0)) //
//     + ( a*( exp2x*(1.0+K) - exp2Kx) / ( exp2Kx*(1.0+exp4x) - 2.0*exp2x) ) //
//     + ( ( a+b-n)*( exp2x*(1.0+K)+exp2Kx ) / ( exp2Kx*(1.0+exp4x) + 2.0*exp2x) );

//   float exp_2x_2Kx = exp2x*exp2Kx;//powf(E,2*(1 + K)*x);
//   float exp_4x_2Kx = exp4x*exp2Kx;//powf(E,2*(2 + K)*x);4x+4Kx
//   float exp_6x_2Kx = exp4x*exp_2x_2Kx; //powf(E,2*(3 + K)*x);
  
    
//   float derivval =
//   (-4*b*exp4x)/powf(-1 + exp4x,2) + 
//    (2*a*(2*exp2x - exp2Kx)*(1 + K))/
//     (-2*exp2x + exp2Kx + exp_4x_2Kx) - 
//    (a*(-exp_2x_2Kx + exp4x*(1 + K))*
//       (-8*exp4x + 2*exp_2x_2Kx*(1 + K) + 2*exp_6x_2Kx*(3 + K)))/
//     powf(-2*exp4x + exp_2x_2Kx + exp_6x_2Kx,2) + 
//    (2*(2*exp2x + exp2Kx)*(1 + K)*(a + b - n))/
//     (2*exp2x + exp2Kx + exp_4x_2Kx) - 
//    (2*(exp_2x_2Kx + exp4x*(1 + K))*
//       (4*exp4x + exp_2x_2Kx*(1 + K) + exp_6x_2Kx*(3 + K))*(a + b - n))/
//     powf(2*exp4x + exp_2x_2Kx + exp_6x_2Kx,2);


//   return VAL/derivval;
// }

// static float
// SLOW_NR(float x,
//         float K,//fix ratio
//         float a,//observed transitions
//         float b,//observed transversions
//         float n){
//   float E = exp(1);
  
//   float fval =(b/(-1 + pow(E,4*x)) + (a*(-pow(E,2*(1 + K)*x) + pow(E,4*x)*(1 + K)))/
//       (-2*pow(E,4*x) + pow(E,2*(1 + K)*x) + pow(E,2*(3 + K)*x)) + 
//      ((pow(E,2*(1 + K)*x) + pow(E,4*x)*(1 + K))*(a + b - n))/
//                (2*pow(E,4*x) + pow(E,2*(1 + K)*x) + pow(E,2*(3 + K)*x)));

//   float derval= 
//    ((-4*b*pow(E,4*x))/pow(-1 + pow(E,4*x),2) + 
//      (2*a*(2*pow(E,2*x) - pow(E,2*K*x))*(1 + K))/
//       (-2*pow(E,2*x) + pow(E,2*K*x) + pow(E,2*(2 + K)*x)) - 
//      (a*(-pow(E,2*(1 + K)*x) + pow(E,4*x)*(1 + K))*
//         (-8*pow(E,4*x) + 2*pow(E,2*(1 + K)*x)*(1 + K) + 2*pow(E,2*(3 + K)*x)*(3 + K)))/
//       pow(-2*pow(E,4*x) + pow(E,2*(1 + K)*x) + pow(E,2*(3 + K)*x),2) + 
//      (2*(2*pow(E,2*x) + pow(E,2*K*x))*(1 + K)*(a + b - n))/
//       (2*pow(E,2*x) + pow(E,2*K*x) + pow(E,2*(2 + K)*x)) - 
//      (2*(pow(E,2*(1 + K)*x) + pow(E,4*x)*(1 + K))*
//         (4*pow(E,4*x) + pow(E,2*(1 + K)*x)*(1 + K) + pow(E,2*(3 + K)*x)*(3 + K))*(a + b - n))/
//     pow(2*pow(E,4*x) + pow(E,2*(1 + K)*x) + pow(E,2*(3 + K)*x),2));

//   cout << " fval " << fval << "   derval " << derval << endl;

//   return fval/derval;
// }

// static int maxiter=0;
// static int numabove =0;
// static int total=0;


// //Binary search for the zero of the derivative of the likelihood
// //which is the same as finding the zero of partof_derivative_likelihood.
// //This function is positive from [0->searchedval] and negative [searchedval->...].

// static float
// _binary_search(
//               float K,//fix ratio
//               float a,//observed transitions
//               float b,//observed transversions
//               float n){
//   float bt_prob = (a/K + b/2)/(n);
//   //cout << "initial " << bt_prob << "   ilen " << interval_length << endl;
//   int i =0;

//   float interval_length = bt_prob;
//   while (  partof_derivative_likelihood(bt_prob,K,a,b,n) > 0 ){
//     i++;
//     bt_prob += interval_length;
//     //cout << "initial " << bt_prob << "   ilen " << interval_length << endl;
//     //PENDING LOOP COUNTER!
//     cout << " i " << endl;
//   }
  

//   //Now we know that the search value is between [bt_prob - interval_length, bt_prob]
//   // BEGIN BINARY SEARCH
//   //uint i = 0;
  
//   while ( interval_length > 0.000001 ){
//     //PENDING LOOP COUNTER!
//     interval_length = interval_length*0.5;
//     float mid = bt_prob - interval_length;
//     if ( partof_derivative_likelihood(mid,K,a,b,n) < 0 )
//       bt_prob = mid;

//     //cout << "ITER " <<(i++) << "  "<< bt_prob << "   ilen " << interval_length << "    dist: " << (K+2)*bt_prob << endl;
//     i++;
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

//   return bt_prob;
// }




// static float
// _NEWTON_RAPHSON(
//               float K,//fix ratio
//               float a,//observed transitions
//               float b,//observed transversions
//               float n){
//   //the initial guess for the number of transversions (i.e. B*t) is
//   // b might be very wrong from expected number of observed.
//   float bt_prob = ((a/(K/2) + b))/(4.0*n);
//   //cout << "initial " << bt_prob << "   ilen " << interval_length << endl;
//   int i =0;

//   float old_bt_prob;
//   do{
//     i++;
//     cout << bt_prob<<endl;
//     old_bt_prob = bt_prob;
//     //bt_prob = old_bt_prob-NEWTON_RAPHSON_RATIO(bt_prob,K,a,b,n);
//     bt_prob = old_bt_prob-SLOW_NR(bt_prob,K,a,b,n);
//   }while ( fabs(old_bt_prob - bt_prob) < 0.000001 && i <20);
//   cout << bt_prob<<endl;
//   if ( i == 20 ){
//     cerr << "WARNING NEWTON RAPHSON FAILED more than " << i << " iterations.\n" << endl;
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

//   return bt_prob;
// }
