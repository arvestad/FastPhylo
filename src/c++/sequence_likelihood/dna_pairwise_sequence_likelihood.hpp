//--------------------------------------------------
//                                        
// File: dna_pairwise_sequence_likelihood.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: dna_pairwise_sequence_likelihood.hpp,v 1.4 2006/12/09 11:45:26 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef DNA_PAIRWISE_SEQUENCE_LIKELIHOOD_HPP
#define DNA_PAIRWISE_SEQUENCE_LIKELIHOOD_HPP

//
// This file contains functions for estimating the ML distance under some different models.
//
// 1. Uncorrected p-distance
// 2. Jukes Cantor
// 3. Kimura 2 parameter
// 4. Tamura Nei
// 5. Felsenstein 84 
//
// Each function takes a struct that provides the sufficient statistica for computing the ML estimate of the
// model parameters. 
//

#include <string>
#include <iostream>
#include <math.h>


//MODEL
typedef enum {
  HAMMING_DISTANCE, //regular hamming distance
  P_DISTANCE,
  JC,
  K2P,
  FAKE_F84,
  TN93
} sequence_model;

//------------------------------
//COMPUTE LIKELIHOOD

static __inline double
compute_JC_log_likelihood_max05(double hamdist, double slen){
  double optprob = hamdist/slen;
  optprob = ( 0.499 < optprob ? 0.499 : optprob );//we don't allow for more than 0.499 probability of change
  if ( hamdist != 0.0 ){
    return hamdist*log(optprob/3.0) + (slen - hamdist)*log(1 - optprob);
  }
  return 0.0;
}


//-------------------------------
// SIMPLE DISTANCE
//
// These structs contain all the information that the correction formulas require
// to compute the corrected distance.
struct simple_string_distance{
  //the number of positions that have been dissregarded due to unknown
  //nucleotide
  int deletedPositions;
  float transitions; 
  float transversions;
};


//TAMURA NEI
struct TN_string_distance{
  int deletedPositions;
  float purine_transitions;
  float pyrimidine_transitions;
  float transversions;
};

// ML STRING DISTANCE
// This is the information that is returned from each computation.
struct ML_string_distance{
  float distance;

  //---------------
  //Change observation probabilities
  float A_A; float A_C; float A_G; float A_T;
  float C_C; float C_G; float C_T;
  float G_G; float G_T;
  float T_T;  
  
};

//---------------------------------------
// utility functions for the structs
simple_string_distance
convert_TN_string_distance_to_simple(TN_string_distance &d);

bool
operator==(const simple_string_distance &d1, const simple_string_distance &d2);

bool
operator==(const TN_string_distance &d1, const TN_string_distance &d2);

std::ostream&
operator<<(std::ostream & os, const simple_string_distance &d);

std::ostream&
operator<<(std::ostream & os, const TN_string_distance &d);

std::ostream&
operator<<(std::ostream & os, const ML_string_distance &sd);

//-------------------------------------------------------------



//-------------------------------------------
//HAMMING DISTANCE
static __inline float
compute_Hamming_distance(simple_string_distance sd){
  return sd.transitions+sd.transversions;
}


//-------------------------------------------
//UNCORRECTED DISTANCE aka p-distance
static __inline float
compute_p_distance(int strlen, simple_string_distance sd){
  // 1 - (len-H)/len i.e. Hamming distance divided by the string length
  return compute_Hamming_distance(sd)/(strlen-sd.deletedPositions);
}

//p-distance that takes the number of gaps into account
static __inline float
compute_p_distance_GAPS(int strlen, simple_string_distance sd, int gaps, float gap_weight){
  // 1 - (len-H)/len i.e. Hamming distance divided by the string length
  return compute_Hamming_distance(sd)/(strlen+gaps*gap_weight);
}


//------------------------------------------------
// JUKES CANTOR

static __inline ML_string_distance
compute_JC(int strlen, simple_string_distance sd){
  ML_string_distance tp;
  tp.distance = -(3.0/4.0) * log( 1 - (4.0/3.0) * compute_p_distance(strlen,sd));

  strlen = strlen-sd.deletedPositions;
  float idprob = ((float)strlen - sd.transitions - sd.transversions)/strlen;
  float changeprob = (1.0- idprob)/3.0;
  
  tp.A_A = idprob; tp.A_C = changeprob;  tp.A_G = changeprob; tp.A_T = changeprob;
  tp.C_C = idprob;  tp.C_G = changeprob; tp.C_T = changeprob;
  tp.G_G = idprob;  tp.G_T = changeprob;
  tp.T_T = idprob;

  return tp;
}

//JC distance that takes the number of gaps into account
static __inline float
compute_JC_GAPS(int strlen, simple_string_distance sd, int gaps, float gap_weight){
  return -(3.0/4.0) * log( 1 - (4.0/3.0) * compute_p_distance_GAPS(strlen, sd,gaps,gap_weight));
}



//---------------------------------------------------
// KIMURA 2 PARAMETER - regular 
//
static __inline ML_string_distance
compute_K2P(int strlen, simple_string_distance sd){

  strlen = strlen-sd.deletedPositions;
  ML_string_distance tp;
  
  float tsprob = sd.transitions/strlen;
  //std::cout << "P... " << sd.transitions <<"/" << strlen << " .... " << P << std::endl;
  float tvprob = sd.transversions/strlen;
  //std::cout << "Q... " << sd.transversions <<"/" << strlen << " .... " << Q << std::endl;

  tp.distance = - 0.5 * log( (1.0 - 2.0*tsprob-tvprob)*sqrt(1.0-2.0*tvprob) );
  //the same as
  //return 0.5*log( 1.0 /(1.0-2.0*P-Q) ) + 1.0/4.0*log( 1.0 /(1.0- 2.0*Q) );
  
  float idprob = ((float)strlen - sd.transitions - sd.transversions)/strlen;
  
  tp.A_A = idprob; tp.A_C = tvprob;  tp.A_G = tsprob; tp.A_T = tvprob;
  tp.C_C = idprob;  tp.C_G = tvprob; tp.C_T = tsprob;
  tp.G_G = idprob;  tp.G_T = tvprob;
  tp.T_T = idprob;


  return tp;
}

// KIMURA 2 PARAMETER WITH FIXED TRANSITION TRANSVERSION RATIO
//
// fixRatio is the A/2B where A is the Markov-transition probability
// of transition and B is the markov-transition probability of
// transversion.  The 2 is there since there are twice as many
// transversions as there are transitions.

ML_string_distance
compute_K2P_fixratio(int strlen, simple_string_distance sd, float fixRatio);


//-----------------------------------------------
//TAMURA NEI - regular
static __inline ML_string_distance
compute_Tamura_Nei(int strlen, TN_string_distance sd,
                   int numAs, int numCs,
                   int numGs, int numTs){

  strlen = strlen - sd.deletedPositions;
  
  float norm = numAs + numCs + numGs + numTs;
  
  float piA = 0.000001+((float)numAs)/norm;
  float piC = 0.000001+((float)numCs)/norm;
  float piG = 0.000001+((float)numGs)/norm;
  float piT = 0.000001+((float)numTs)/norm;


  float piR = piA+piG;
  float piY = piC+piT;


  ML_string_distance tp;
  
  float ts_purine_prob = sd.purine_transitions/strlen;
  float ts_pyrimidine_prob = sd.pyrimidine_transitions/strlen;
  float tv_prob = sd.transversions/strlen;
  
  tp.distance = -2.0*piA*piG/piR*log(1.0-piR*ts_purine_prob/(2.0*piA*piG) - tv_prob/(2.0*piR))
     -2.0*piT*piC/piY*log(1.0-piY*ts_pyrimidine_prob/(2.0*piC*piT) - tv_prob/(2.0*piY))
    - 2.0*( piR*piY - piA*piG*piY/piR - piT*piC*piR/piY)*log(1.0-tv_prob/(2.0*piR*piY));

  
  tp.A_A = 1.0-ts_purine_prob-tv_prob; tp.A_C = tv_prob*0.5;  tp.A_G = ts_purine_prob; tp.A_T = tv_prob*0.5;
  tp.C_C = 1.0-ts_pyrimidine_prob-tv_prob;  tp.C_G = tv_prob*0.5; tp.C_T = ts_pyrimidine_prob;
  tp.G_G = 1.0-ts_purine_prob-tv_prob;  tp.G_T = tv_prob*0.5;
  tp.T_T = 1.0-ts_pyrimidine_prob-tv_prob;


  return tp;

}

// TAMURA NEI WITH FIX RATIO
ML_string_distance
compute_Tamura_Nei_fixratio(int strlen, TN_string_distance sd,
                            int numAs, int numCs,
                            int numGs, int numTs,
                            float purine_ts_tv_ratio,
                            float pyrimidine_ts_tv_ratio);



//---------------------------------------------------
//FAKE F84
//
//WARNING
//These methods proved some kind of estimate for the F84 model but it
//is not the ML estimate.
typedef struct{
  float A;
  float B;
  float C;
} FAKE_F84_constants;

static __inline FAKE_F84_constants
compute_FAKE_F84_constants(int basefreqA,
                      int basefreqC,
                      int basefreqG,
                      int basefreqT){
  FAKE_F84_constants consts;
  int len = basefreqA + basefreqC + basefreqG + basefreqT;
  float piA = 0.000001+((float)basefreqA)/len;
  float piC = 0.000001+((float)basefreqC)/len;
  float piG = 0.000001+((float)basefreqG)/len;
  float piT = 0.000001+((float)basefreqT)/len;


  float piY = piC+piT;
  float piR = piA+piG;
  consts.A = piC*piT/piY + piA*piG/piR;
  consts.B = piC*piT + piA*piG;
  consts.C = piR*piY;
  return consts;
}


static __inline ML_string_distance
compute_FAKE_F84(int strlen, simple_string_distance sd,
            FAKE_F84_constants consts){

  ML_string_distance tp;
  strlen = strlen -sd.deletedPositions;

  float P = sd.transitions/strlen;
  float Q = sd.transversions/strlen;

  
  //PENDING REALLY SLOW
  tp.distance = -2.0*consts.A*log( 1.0  - P / (2.0*consts.A) - (consts.A-consts.B)*Q/(2.0*consts.A*consts.C) )
    + 2.0*(consts.A - consts.B - consts.C)*log(1.0 - Q/(2.0*consts.C));

   float idprob = ((float)strlen - sd.transitions - sd.transversions)/strlen;
  float tsprob = ((float) sd.transversions)/strlen;
  float tvprob = (1.0 - idprob -tsprob)/2.0;
  
  tp.A_A = idprob; tp.A_C = tvprob;  tp.A_G = tsprob; tp.A_T = tvprob;
  tp.C_C = idprob;  tp.C_G = tvprob; tp.C_T = tsprob;
  tp.G_G = idprob;  tp.G_T = tvprob;
  tp.T_T = idprob;

  return tp;
}





#endif // DNA_PAIRWISE_SEQUENCE_LIKELIHOOD_HPP
