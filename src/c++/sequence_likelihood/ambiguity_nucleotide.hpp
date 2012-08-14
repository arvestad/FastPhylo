//--------------------------------------------------
//                                        
// File: ambiguity_nucleotide.hpp
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: ambiguity_nucleotide.hpp,v 1.21 2006/12/29 09:07:10 isaac Exp $
//
//--------------------------------------------------

#ifndef AMBIGUITY_NUCLEOTIDE_HPP
#define AMBIGUITY_NUCLEOTIDE_HPP

#include "nucleotide.hpp"
#include "dna_pairwise_sequence_likelihood.hpp"
#include <assert.h>
#include <string>
#include <cstdlib>
#include <cstdlib>
#include <iostream>
#include "log_utils.hpp"

//AMBIGUITY_NUCLEOTIDE
typedef struct{
  nucleotide n;
  //probabilties for each nucleotide
  float probA;
  float probC;
  float probG;
  float probT;
} ambiguity_nucleotide;

//DISTANCE
typedef struct{
  float purine_transition_prob;
  float pyrimidine_transition_prob;
  float transversion_prob;
} ambiguity_distance;


std::ostream&
operator<<(std::ostream & os, const ambiguity_distance dist);

static __inline bool
operator==(const ambiguity_nucleotide &a1, const ambiguity_nucleotide &a2){

  return ( a1.n==a2.n && a1.probA == a2.probA && a1.probC == a2.probC
           && a1.probG == a2.probG && a1.probT == a2.probT );
}

static __inline bool
operator!=(const ambiguity_nucleotide &a1, const ambiguity_nucleotide &a2){
  return !(a1==a2);
}
//----------------------------------
//CONVERSION nucleotide -> ambiguity_nucleotide
//
// This function assumes uniform distribution of the nucleotides.  To
// avoid overflow of probabilities probability 0.999999 is used
// instead of 1.0 at some points.
//
// WARNING! DNA_UNKOWN_ and DNA_NOT_ALLOWED are not allowed inputs to these
// functions.
static __inline ambiguity_nucleotide
nucleotide2ambiguity_nucleotideUNIFORM(const nucleotide &n){
  
  switch (n){//                                  A     C    G    T
  case DNA_A_: { return (ambiguity_nucleotide){n, 1.0,   0,   0,   0}; }
  case DNA_C_: { return (ambiguity_nucleotide){n,   0, 1.0,   0,   0}; }
  case DNA_G_: { return (ambiguity_nucleotide){n,   0,   0, 1.0,   0}; }
  case DNA_T_: { return (ambiguity_nucleotide){n,   0,   0,   0, 1.0}; }
  case DNA_M_: { return (ambiguity_nucleotide){n, 0.5, 0.5,   0,   0}; }
  case DNA_R_: { return (ambiguity_nucleotide){n, 0.5,   0, 0.5,   0}; }//an.probA = 0.5; an.probG = 0.5; break;
  case DNA_W_: { return (ambiguity_nucleotide){n, 0.5,   0,   0, 0.5}; }//an.probA = 0.5; an.probT = 0.5; break;
  case DNA_S_: { return (ambiguity_nucleotide){n,   0, 0.5, 0.5,   0}; }//an.probC = 0.5; an.probG = 0.5; break;
  case DNA_Y_: { return (ambiguity_nucleotide){n,   0, 0.5,   0, 0.5}; }//an.probC = 0.5; an.probT = 0.5; break;
  case DNA_K_: { return (ambiguity_nucleotide){n,   0,   0, 0.5, 0.5}; }//an.probG = 0.5; an.probT = 0.5; break;
  case DNA_V_: { return (ambiguity_nucleotide){n,0.333333,0.333333,0.333333, 0}; }//an.probA = 0.999999/3.0; an.probC = 0.999999/3.0; an.probG = 0.999999/3.0; break;
  case DNA_H_: { return (ambiguity_nucleotide){n, 0.333333, 0.333333, 0, 0.333333}; }//an.probA = 0.999999/3.0; an.probC = 0.999999/3.0; an.probT = 0.999999/3.0; break;
  case DNA_D_: { return (ambiguity_nucleotide){n, 0.333333, 0, 0.333333, 0.333333}; }//an.probA = 0.999999/3.0; an.probG = 0.999999/3.0; an.probT = 0.999999/3.0; break;
  case DNA_B_: { return (ambiguity_nucleotide){n, 0, 0.333333, 0.333333, 0.333333}; }//an.probC = 0.999999/3.0; an.probG = 0.999999/3.0; an.probT = 0.999999/3.0; break;
  case DNA_N_: { return (ambiguity_nucleotide){n, 0.25, 0.25, 0.25, 0.25}; }//an.probA = 0.25; an.probC = 0.25; an.probG = 0.25; an.probT = 0.25; break;
  case DNA_NOT_ALLOWED: USER_WARNING( "DNA_NOT_ALLOWED can not be used as input" ); break;
  case DNA_UNKNOWN_: USER_WARNING( "DNA_UNKOWN can not be used as input" ); break; //unkown is not handled
    //  default: USER_WARNING( "AMBIG_X_FLAG can not be used as input: "<<n );
  }
  
   return (ambiguity_nucleotide){n,0,0,0,0};
}

//assuming that the nucleotid is a regular
static __inline ambiguity_nucleotide
regularnucleotide2ambiguity_nucleotide(const nucleotide &n){
  
  switch (n){
  case DNA_A_: { return (ambiguity_nucleotide){n, 1.0,   0,   0,   0}; }
  case DNA_C_: { return (ambiguity_nucleotide){n,   0, 1.0,   0,   0}; }
  case DNA_G_: { return (ambiguity_nucleotide){n,   0,   0, 1.0,   0}; }
  case DNA_T_: { return (ambiguity_nucleotide){n,   0,   0,   0, 1.0}; }
  default:
    USER_WARNING("Input has to be regular nucleotide"); break;
  }

  
  return (ambiguity_nucleotide){n,0,0,0,0};
}


  // CONVERSION WITH BASE FREQUENCES
static __inline ambiguity_nucleotide
nucleotide2ambiguity_nucleotide(const nucleotide &n,
                                int basefreqA = 1,
                                int basefreqC = 1,
                                int basefreqG = 1,
                                int basefreqT = 1){
  ambiguity_nucleotide an = {n,0,0,0,0};
  
  switch (n){
  case DNA_A_: an.probA = 1.0; break;
  case DNA_C_: an.probC = 1.0; break;
  case DNA_G_: an.probG = 1.0; break;
  case DNA_T_: an.probT = 1.0; break;
  case DNA_UNKNOWN_: USER_ERROR("DNA_UNKNOWN_ can not be input"); break; //unkown is not handled
  case DNA_M_: an.probA = (0.999999*basefreqA)/(basefreqA+basefreqC); an.probC = (0.999999*basefreqC)/(basefreqA+basefreqC); break;
  case DNA_R_: an.probA = (0.999999*basefreqA)/(basefreqA+basefreqG); an.probG = (0.999999*basefreqG)/(basefreqA+basefreqG); break;
  case DNA_W_: an.probA = (0.999999*basefreqA)/(basefreqA+basefreqT); an.probT = (0.999999*basefreqT)/(basefreqA+basefreqT); break;
  case DNA_S_: an.probC = (0.999999*basefreqC)/(basefreqC+basefreqG); an.probG = (0.999999*basefreqG)/(basefreqC+basefreqG); break;
  case DNA_Y_: an.probC = (0.999999*basefreqC)/(basefreqC+basefreqT); an.probT = (0.999999*basefreqT)/(basefreqC+basefreqT); break;
  case DNA_K_: an.probG = (0.999999*basefreqG)/(basefreqG+basefreqT); an.probT = (0.999999*basefreqT)/(basefreqG+basefreqT); break;
  case DNA_V_:
    an.probA = (0.999999*basefreqA)/(basefreqA+basefreqC+basefreqG);
    an.probC = (0.999999*basefreqC)/(basefreqA+basefreqC+basefreqG);
    an.probG = (0.999999*basefreqG)/(basefreqA+basefreqC+basefreqG);
    break;
  case DNA_H_:
    an.probA = (0.999999*basefreqA)/(basefreqA+basefreqC+basefreqT);
    an.probC = (0.999999*basefreqC)/(basefreqA+basefreqC+basefreqT);
    an.probT = (0.999999*basefreqT)/(basefreqA+basefreqC+basefreqT);
    break;
  case DNA_D_:
    an.probA = (0.999999*basefreqA)/(basefreqA+basefreqG+basefreqT);
    an.probG = (0.999999*basefreqG)/(basefreqA+basefreqG+basefreqT);
    an.probT = (0.999999*basefreqT)/(basefreqA+basefreqG+basefreqT);
    break;
  case DNA_B_:
    an.probC = (0.999999*basefreqC)/(basefreqC+basefreqG+basefreqT);
    an.probG = (0.999999*basefreqG)/(basefreqC+basefreqG+basefreqT);
    an.probT = (0.999999*basefreqT)/(basefreqC+basefreqG+basefreqT);
    break;
  case DNA_N_:
    an.probA = (0.999999*basefreqA)/(basefreqA+basefreqC+basefreqG+basefreqT);
    an.probC = (0.999999*basefreqC)/(basefreqA+basefreqC+basefreqG+basefreqT);
    an.probG = (0.999999*basefreqG)/(basefreqA+basefreqC+basefreqG+basefreqT);
    an.probT = (0.999999*basefreqT)/(basefreqA+basefreqC+basefreqG+basefreqT);
    break;
  case DNA_NOT_ALLOWED: USER_WARNING("DNA_NOT_ALLOWED can not be input" ); break;
  default: USER_WARNING("AMBIG_X_FLAG can not be input: " << n); 
  }
  return an;
}


//--------------------------------------------
// AMBIGUITY DISTANCE
//
// Computes the probability of transition and transversion between the
// two nucleotides.
static __inline ambiguity_distance
compute_ambiguity_distance(const ambiguity_nucleotide &a1, const ambiguity_nucleotide &a2){

  float matchprob = a1.probA * a2.probA;
  matchprob += a1.probC * a2.probC;
  matchprob += a1.probG * a2.probG;
  matchprob += a1.probT * a2.probT;
  
  ambiguity_distance dist = { a1.probA*a2.probG + a1.probG*a2.probA,//purine trans
                              a1.probC*a2.probT + a1.probT*a2.probC,//pyrimidine trans
                              0};//tv
  dist.transversion_prob = 1.0-matchprob-dist.purine_transition_prob-dist.pyrimidine_transition_prob;

  //PENDING slow with all these tests
  if ( dist.transversion_prob < 0 )
    dist.transversion_prob = 0;
  if ( dist.purine_transition_prob < 0 )
    dist.purine_transition_prob = 0;
  if (dist.pyrimidine_transition_prob < 0 )
    dist.pyrimidine_transition_prob = 0;
  if ( matchprob < 0 )
    matchprob = 0;
  float total = matchprob + dist.purine_transition_prob +dist.pyrimidine_transition_prob  + dist.transversion_prob;
  matchprob = matchprob/total;
  dist.transversion_prob = dist.transversion_prob/total;
  dist.purine_transition_prob = dist.purine_transition_prob/total;
  dist.pyrimidine_transition_prob = dist.pyrimidine_transition_prob/total;

  //#warning "ambig if statements"
  // //things that need to be true
//   assert ( matchprob <= 1.0 );
//   assert ( matchprob >= 0.0 );
//   assert ( dist.purine_transition_prob <= 1.0 );
//   assert ( dist.purine_transition_prob >= 0.0 );
//   assert ( dist.pyrimidine_transition_prob <= 1.0 );
//   assert ( dist.pyrimidine_transition_prob >= 0.0 );
//   assert ( dist.transversion_prob <= 1.0 );
//   assert ( dist.transversion_prob >= 0.0 );
//   assert ( 0.999995 <= matchprob + dist.purine_transition_prob +dist.pyrimidine_transition_prob + dist.transversion_prob );
//   assert ( 1.0 >= matchprob + dist.purine_transition_prob +dist.pyrimidine_transition_prob  + dist.transversion_prob);
  
  return dist;
}

//using transition probabilities
static __inline ambiguity_distance
compute_ambiguity_distance_using_transition_probabilities(const ambiguity_nucleotide &a1, const ambiguity_nucleotide &a2,
							  const ML_string_distance &tp){

//   float matchprob = a1.probA * a2.probA*tp.A_A
//     +a1.probC * a2.probC*tp.C_C
//     +a1.probG * a2.probG*tp.G_G
//     + a1.probT * a2.probT*tp.T_T;
//   //std::cout << "matchprob " << matchprob <<std::endl;
//   ambiguity_distance dist = { (a1.probA*a2.probG + a1.probG*a2.probA)*tp.A_G,//purine
//                               (a1.probC*a2.probT + a1.probT*a2.probC)*tp.C_T,//pyrimidine
//                               0};//tv
//   //std::cout << "purine " << dist.purine_transition_prob <<std::endl;
//   //std::cout << "pyrimidine " << dist.pyrimidine_transition_prob <<std::endl;
  
//   dist.transversion_prob = a1.probA*a2.probC*tp.A_C + a1.probA*a2.probT*tp.A_T +
// 			    a1.probC*a2.probA*tp.A_C + a1.probT*a2.probA*tp.A_T +
// 			    a1.probG*a2.probC*tp.C_G + a1.probG*a2.probT*tp.G_T +
// 			    a1.probC*a2.probG*tp.C_G + a1.probT*a2.probG*tp.G_T;
//   //std::cout << "transv " << dist.transversion_prob <<std::endl;
  //SPEED UP
  float matchprob;
  ambiguity_distance dist;

  float a1prob = a1.probA;
  const float a2probA = a2.probA;
  const float a2probC = a2.probC;
  const float a2probG = a2.probG;
  const float a2probT = a2.probT;
  
  // A->
  matchprob = a1prob*a2probA*tp.A_A;
  dist.purine_transition_prob = a1prob * a2probG*tp.A_G;
  dist.transversion_prob = a1prob * (a2probC*tp.A_C + a2probT*tp.A_T);
  
  // C->
  a1prob = a1.probC;
  matchprob += a1prob*a2probC*tp.C_C;
  dist.pyrimidine_transition_prob = a1prob * a2probT*tp.C_T;
  dist.transversion_prob += a1prob * (a2probA*tp.A_C + a2probG*tp.C_G);

  // G->
  a1prob = a1.probG;
  matchprob += a1prob*a2probG*tp.G_G;
  dist.purine_transition_prob += a1prob * a2probA*tp.A_G;
  dist.transversion_prob += a1prob * (a2probC*tp.C_G + a2probT*tp.G_T);
  
  // T->
  a1prob = a1.probT;
  matchprob += a1prob*a2probT*tp.T_T;
  dist.pyrimidine_transition_prob += a1prob * a2probC*tp.C_T;
  dist.transversion_prob += a1prob * (a2probA*tp.A_T + a2probG*tp.G_T);
  
  //-----------------------------------
  float totalprob = matchprob + 
    dist.purine_transition_prob + dist.pyrimidine_transition_prob + dist.transversion_prob;

  //if totalprob is 0 this means that the transition probabilities are zero for every
  //possible transition between the ambigious symbols. In particular this means that the ambious symbols
  //dont match. So we recompute the transition and transversions without the transition probs.
  if ( totalprob == 0 ){
    dist.purine_transition_prob = (a1.probA*a2.probG + a1.probG*a2.probA);//purine
    dist.pyrimidine_transition_prob = (a1.probC*a2.probT + a1.probT*a2.probC);//pyrimidine
                                
    dist.transversion_prob = a1.probA*a2.probC + a1.probA*a2.probT +
    a1.probC*a2.probA + a1.probT*a2.probA +
    a1.probG*a2.probC + a1.probG*a2.probT +
     a1.probC*a2.probG + a1.probT*a2.probG;
    matchprob = 0;
    totalprob = matchprob + dist.purine_transition_prob + dist.pyrimidine_transition_prob + dist.transversion_prob;
  }

  dist.purine_transition_prob = dist.purine_transition_prob/totalprob;
  dist.pyrimidine_transition_prob = dist.pyrimidine_transition_prob/totalprob;
  dist.transversion_prob = dist.transversion_prob/totalprob;
  matchprob = matchprob/totalprob;//pending unnecesary
  
  //std::cout << "sum " << matchprob + dist.purine_transition_prob +dist.pyrimidine_transition_prob + dist.transversion_prob << std::endl;
  //things that need to be true
  assert ( matchprob <= 1.0 );
  assert ( matchprob >= 0.0 );
  assert ( dist.purine_transition_prob <= 1.0 );
  assert ( dist.purine_transition_prob >= 0.0 );
  assert ( dist.pyrimidine_transition_prob <= 1.0 );
  assert ( dist.pyrimidine_transition_prob >= 0.0 );
  assert ( dist.transversion_prob <= 1.0 );
  assert ( dist.transversion_prob >= 0.0 );
  assert ( 0.999995 < matchprob + dist.purine_transition_prob +dist.pyrimidine_transition_prob + dist.transversion_prob );
  assert ( 1.000001 >= matchprob + dist.purine_transition_prob +dist.pyrimidine_transition_prob  + dist.transversion_prob);
  
  return dist;
}

//
// Does the same as above but it can only take a regular nucleotide A,G,C,T
//
static __inline ambiguity_distance
compute_ambiguity_distance_using_transition_probabilities(const nucleotide &n, const ambiguity_nucleotide &a2,
							  const ML_string_distance &tp){
  //SPEED UP
  float matchprob;
  ambiguity_distance dist = {0,0,0};

  const float a2probA = a2.probA;
  const float a2probC = a2.probC;
  const float a2probG = a2.probG;
  const float a2probT = a2.probT;
  float totalprob=0;

  switch(n){
  case DNA_A_:
    // A->
    matchprob = a2probA*tp.A_A;
    dist.purine_transition_prob = a2probG*tp.A_G;
    dist.transversion_prob = (a2probC*tp.A_C + a2probT*tp.A_T);
    
    totalprob = matchprob + dist.purine_transition_prob + dist.transversion_prob;
    if (totalprob == 0){
      matchprob = a2probA;
      dist.purine_transition_prob = a2probG;
      dist.transversion_prob = (a2probC + a2probT);
      totalprob = matchprob + dist.purine_transition_prob + dist.transversion_prob;
    }
    break;
  case DNA_C_:
    // C->
    matchprob = a2probC*tp.C_C;
    dist.pyrimidine_transition_prob =  a2probT*tp.C_T;
    dist.transversion_prob = (a2probA*tp.A_C + a2probG*tp.C_G);
    
    totalprob = matchprob + dist.pyrimidine_transition_prob + dist.transversion_prob;
    if (totalprob == 0){
      matchprob = a2probC;
      dist.pyrimidine_transition_prob =  a2probT;
      dist.transversion_prob = (a2probA + a2probG);
      totalprob = matchprob + dist.pyrimidine_transition_prob + dist.transversion_prob;
    }
    break;
  case DNA_G_:
    // G->
    matchprob = a2probG*tp.G_G;
    dist.purine_transition_prob =  a2probA*tp.A_G;
    dist.transversion_prob = (a2probC*tp.C_G + a2probT*tp.G_T);
    
    totalprob = matchprob + dist.purine_transition_prob + dist.transversion_prob;
    if (totalprob == 0){
      matchprob = a2probG;
      dist.purine_transition_prob =  a2probA;
      dist.transversion_prob = (a2probC + a2probT);
      totalprob = matchprob + dist.purine_transition_prob + dist.transversion_prob;
    }
    break;
  case DNA_T_:
    // T->
    matchprob = a2probT*tp.T_T;
    dist.pyrimidine_transition_prob =  a2probC*tp.C_T;
    dist.transversion_prob = (a2probA*tp.A_T + a2probG*tp.G_T);
    
    totalprob = matchprob + dist.pyrimidine_transition_prob + dist.transversion_prob;
    if (totalprob == 0){
      matchprob = a2probT;
      dist.pyrimidine_transition_prob =  a2probC;
      dist.transversion_prob = (a2probA + a2probG);
      totalprob = matchprob + dist.pyrimidine_transition_prob + dist.transversion_prob;
    }
    break;
  default:
    USER_WARNING("Shouldn't be here: can only use regular nucleotides");
    return dist;
  }
  //-----------------------------------
  
  dist.purine_transition_prob = dist.purine_transition_prob/totalprob;
  dist.pyrimidine_transition_prob = dist.pyrimidine_transition_prob/totalprob;
  dist.transversion_prob = dist.transversion_prob/totalprob;
  
  return dist;
}



//--------------------------------------------------
// Checks if all possible nucleotides of n are possible nucleotides in
// the nucleotide set. E.g. A is in {A,G}, {A,C} is in {A,C,G}, T is in {T}.
//

static __inline bool
is_ambiguity_contained(const ambiguity_nucleotide &n, const ambiguity_nucleotide &set){

  if ( n.probA > 0 && set.probA == 0 )
    return false;
  if ( n.probC > 0 && set.probC == 0 )
    return false;
  if ( n.probG > 0 && set.probG == 0 )
    return false;
  if ( n.probT > 0 && set.probT == 0 )
    return false;
  return true;
}


bool
is_ambiguity_contained(const nucleotide &n_, const nucleotide &set_);

#endif // AMBIGUITY_NUCLEOTIDE_H
