//--------------------------------------------------
//                                        
// File: AML_star.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: AML_star.cpp,v 1.1 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------

#include "AML_star.hpp"
#include "string_compare.hpp"
#include "LikelihoodMatrix.hpp"
#include <string>
#include "log_utils.hpp"
#include "SequenceTree.hpp"
#include "dna_pairwise_sequence_likelihood.hpp"

static char
parsimony_star(char a, char b, char c){
  //if two are the same then return that one.
  if ( a==b )
    return a;

  if ( a==c )
    return a;

  if ( b==c )
    return b;

  //if none equal we always reutrn a (this is used below in the algorithm)
  return a;
}



double 
find_AML_star(std::string &a, 
	 std::string &b, 
	 std::string &c,
	 std::string &center,
         sequence_model model){
  std::vector<std::string *> vec(4);
  
  vec[0] = &a;
  vec[1] = &b;
  vec[2] = &c;
  //add parsimony solution
  center.reserve(a.size());
  center.clear();
  for ( size_t i = 0 ; i < a.size() ; i++ )
    center.push_back(parsimony_star(a[i],b[i],c[i]));
  vec[3] = &center;

  SequenceTree::NodeMatrix lm(4);
  fillLikelihoodMatrix(vec,lm,model);
  
  //likelihood with x in center
  double a_l = lm.getDistance(0,1) + lm.getDistance(0,2);
  double b_l = lm.getDistance(0,1) + lm.getDistance(1,2);
  double c_l = lm.getDistance(0,2) + lm.getDistance(1,2);
  double parsimony_l = lm.getDistance(0,3) + lm.getDistance(1,3) + lm.getDistance(2,3);  
  //  PRINT(a_l); PRINT(b_l); PRINT(c_l); PRINT(parsimony_l);

  double max_l = parsimony_l;
  std::string *maxstr = &center;
  if ( a_l > max_l ){
    //    PRINT(a_l);
    maxstr = &a;
    max_l = a_l;
  }

  if ( b_l > max_l ){
    //PRINT(b_l);
    maxstr = &b;
    max_l = b_l;
  }

  if ( c_l > max_l ){
    //PRINT(c_l);
    maxstr = &c;
    max_l = c_l;
  }

  
  //if parsimony not min copy the correct one
  if ( maxstr != &center ){
    center.clear();
    center.append(*maxstr);
  }

  return max_l;
}

double
find_most_likely_parsimonious(std::string &a, 
                              std::string &b, 
                              std::string &c,
                              std::string &center){

  int numAllDifferent =0;
  int numAparsi = 0;
  int numBparsi = 0;
  int numCparsi = 0;

  int len = a.length();
  for ( int i = 0 ; i < len ; i++ ){
    if ( a[i] == b[i] ){
      numAparsi++;numBparsi++;
      if ( a[i] == c[i] )
        numCparsi++;
    }
    else if ( a[i] == c[i] ){
      numAparsi++;numCparsi++;
    }
    else if ( b[i] == c[i] ){
      numBparsi++;numCparsi++;
    }
    else
      numAllDifferent++;
  }

  //find closest
  std::string *closest = &a;
  if ( numBparsi > numAparsi )
    closest = &b;
  if ( numCparsi > numAparsi )
    closest = &c;

  for ( int i = 0 ; i < len ; i++ ){
    if ( a[i] == b[i] )
      center[i] = a[i];
    else if ( a[i] == c[i] )
      center[i] = a[i];
    else if ( b[i] == c[i] )
      center[i] = b[i];
    else//if all different set to a
      center[i] = (*closest)[i];
  }

  double  likelihood =
    compute_JC_log_likelihood_max05((len-numAparsi-numAllDifferent),len) //a->parsimony
    +compute_JC_log_likelihood_max05((len-numBparsi),len) //b->parsimony
    +compute_JC_log_likelihood_max05((len-numCparsi),len); //c->parsimony

  return likelihood;
}



bool
improve_AML_star(std::string &a, 
                 std::string &b, 
                 std::string &c,
                 std::string &center,
                 sequence_model model){

  if ( model != JC && model != P_DISTANCE )
    PROG_ERROR("not supported");

  
  double slen = 1.0*a.length();


  double hamdist;
  
  //----------------
  //current likelihood
  double currL = 0.0;
  
  //likelihood for a-center
  hamdist = hamming_distance(a,center);
  currL += compute_JC_log_likelihood_max05(hamdist,slen);

  //likelihood for b-center
  hamdist = hamming_distance(b,center);
  currL += compute_JC_log_likelihood_max05(hamdist,slen);
  
  //likelihood for c-center
  hamdist = hamming_distance(c,center);
  currL += compute_JC_log_likelihood_max05(hamdist,slen);
  
  //-----------------
  //a likelihoood
  double aL = 0.0;
  double bL = 0.0;
  double cL = 0.0;
  double tmp;
  
  //a-b
  hamdist = hamming_distance(a,b);
  tmp = compute_JC_log_likelihood_max05(hamdist,slen);
  aL += tmp;  
  bL += tmp;
  
  //a-c
  hamdist = hamming_distance(a,c);
  tmp = compute_JC_log_likelihood_max05(hamdist,slen);
  aL+=tmp;  
  cL+=tmp;
  
  //b-c
  hamdist = hamming_distance(b,c);
  tmp = compute_JC_log_likelihood_max05(hamdist,slen);
  cL += tmp;
  bL += tmp;


  //--------------
  //parsimony
  double parsimonyL=find_most_likely_parsimonious(a,b,c,center);
  
  
  // PRINT(currL);PRINT(aL);PRINT(bL);PRINT(cL);PRINT(parsimonyL);
  //---------------
  //find maximum
  double max_l = parsimonyL;
  std::string *maxstr = &center;
  if ( aL > max_l ){
    maxstr = &a;
    max_l = aL;
  }
  if ( bL > max_l ){
    maxstr = &b;
    max_l = bL;
  }
  if ( cL > max_l ){
    maxstr = &c;
    max_l = cL;
  }

  //if parsimony not min copy the correct one
  if ( maxstr != &center ){
    center.clear();
    center.append(*maxstr);
  }

  if ( max_l - currL > 0.00001 ){
    //PRINT(currL - max_l);PRINT(currL -aL);
    return true;
  }

  return false;
}










