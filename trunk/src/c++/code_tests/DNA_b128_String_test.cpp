//--------------------------------------------------
//                                        
// File: DNAtest.cpp
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: DNA_b128_String_test.cpp,v 1.19 2006/12/29 09:07:10 isaac Exp $                                 
//
//--------------------------------------------------


#include "DNA_b128_String.hpp"
#include "string_compare.hpp"
#include <string>
#include <iostream>
#include <math.h>
#include "log_utils.hpp"

//we want the assertions to work here
#undef NDEBUG

using namespace std;


#define TEST_STRING_DISTANCE(dna1,dna2,str1,str2)\
{\
  PRINT_TIME(d1 = DNA_b128_String::computeDistance(dna1,dna2););\
  PRINT_TIME(t1 = DNA_b128_String::computeTAMURANEIDistance(dna1,dna2););\
  PRINT_TIME(t2 = TN_string_compare(str1,str2););\
  d2 = convert_TN_string_distance_to_simple(t2);\
  cout << d1 << endl;\
  cout << d2 << endl;\
  cout << t1 << endl;\
  cout << t2 << endl;\
  assert (d1 == d2 );\
  assert (t1 == t2 );\
  assert ( str2 == dna2.toString());\
 assert ( str1 == dna1.toString() );\
}


void operator+=(TN_string_distance &sum, const TN_string_distance &add){
  sum.deletedPositions += add.deletedPositions;
  sum.purine_transitions += add.purine_transitions;
  sum.pyrimidine_transitions += add.pyrimidine_transitions;
  sum.transversions += add.transversions;
}

void operator+=(simple_string_distance &sum, const simple_string_distance &add){
  sum.deletedPositions += add.deletedPositions;
  sum.transitions += add.transitions;
  sum.transversions += add.transversions;
}

void
set_empty(simple_string_distance &d){
  d.deletedPositions = 0;
  d.transversions = 0;
  d.transitions = 0;
}

void
set_empty(TN_string_distance &d){
  d.deletedPositions = 0;
  d.purine_transitions = 0;
  d.pyrimidine_transitions =0;
  d.transversions = 0;
}

void
timing_test(size_t times=1000){
  cout <<"--------------\n TIMING TEST   times=" << times<< endl;

  string str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg";//one b128
  str1 = str1+str1+str1+str1;// four b128
  str1 = str1+str1+str1+str1;//i.e. sixteen b128
  str1 = str1+str1+str1+str1;//64 b128s
  str1 = str1+str1+str1+str1;//256 b128s
  str1 = str1+str1+str1+str1;//1024 b128s
  str1 = str1+str1+str1+str1;//4096 b128s
  str1 = str1+str1+str1+str1;//16384 b128s
  str1 = str1+str1+str1+str1;//... b128s
  
  cout <<"creating strings of length " << str1.length() << endl;
  DNA_b128_String dna1(str1.length(),str1);
  string str2 = str1;
  //A->G is a transition A->C a transversion
  for ( unsigned int i = 0 ; i < str1.length()/2 ; i++ ){
    int rand_index = (int) (((float)str2.length())*rand()/(RAND_MAX+1.0));
    str2[rand_index] = 'a';
    if ( i % 5 == 0 )
      str2[rand_index] ='-';
  }

  DNA_b128_String dna2(str2.length(),str2);
  //full compute
  TN_string_distance sum_full_TN;
  set_empty(sum_full_TN);
  clock_t time_full = clock();
  for(size_t time=0; time<times; time++){
    float divergence_matrix[4][4];
    int dels = complete_dna_string_compare(divergence_matrix,str1,str2);
    sum_full_TN += divergence_matrix_2_TN_distance(divergence_matrix,dels);
  }
  time_full = clock() - time_full;
  
  //REGULAR
  TN_string_distance sum_regular_TN;
  simple_string_distance sum_regular_simple;
  set_empty(sum_regular_TN);
  set_empty(sum_regular_simple);
  clock_t time_regular = clock();
  for(size_t time=0; time<times; time++){
    TN_string_distance d = TN_string_compare(str1,str2);
    sum_regular_TN += d;
    sum_regular_simple += convert_TN_string_distance_to_simple(d);
  }
  time_regular = clock() - time_regular;

  //TN compute
  TN_string_distance sum_TN;
  set_empty(sum_TN);
  clock_t time_TN = clock();
  for(size_t time=0; time<times; time++){
    TN_string_distance d = DNA_b128_String::computeTAMURANEIDistance(dna1,dna2);
    sum_TN += d;		         
  }
  time_TN = clock() - time_TN;


  //K2P compute
  simple_string_distance sum_K2P;
  set_empty(sum_K2P);
  clock_t time_K2P = clock();
  for(size_t time=0; time<times; time++){
    simple_string_distance d = DNA_b128_String::computeDistance(dna1,dna2);
    sum_K2P += d;		         
  }
  time_K2P = clock() - time_K2P;

  cout << "FULL TIME       : " << (double(time_full)/(CLOCKS_PER_SEC/1000)) << " ms"<< endl;
  cout << "REGULAR TIME    : " << (double(time_regular)/(CLOCKS_PER_SEC/1000)) << " ms"<< endl;
  cout << "TN TIME         : " << (double(time_TN)/(CLOCKS_PER_SEC/1000)) << " ms"<< endl;
  cout << "K2P TIME        : " << (double(time_K2P)/(CLOCKS_PER_SEC/1000)) << " ms"<< endl;
  
  PRINT(sum_TN);
  PRINT(sum_K2P);
  
  ASSERT_EQ(sum_regular_simple,sum_K2P );
  ASSERT_EQ( sum_regular_TN,sum_TN );
  ASSERT_EQ( sum_regular_TN,sum_full_TN );

  cout << "END TIMING TEST"<< endl;
}


// DNA_M_ = 0x300,// 	A or C 	        K
//   DNA_R_ = 0x500,// 	A or G 	        Y
//   DNA_W_ = 0x900,// 	A or T 	        W
//   DNA_S_ = 0x600,// 	C or G 	        S
//   DNA_Y_ = 0xA00,// 	C or T 	        R
//   DNA_K_ = 0xC00,// 	G or T 	        M
//   DNA_V_ = 0x700,// 	A or C or G 	B
//   DNA_H_ = 0xB00,// 	A or C or T 	D
//   DNA_D_ = 0xD00,// 	A or G or T 	H
//   DNA_B_ = 0xE00,// 	C or G or T 	V
//   DNA_N_ = 0xF00,//  A or C or G or T  // 'N', 'X', 'x', or '?',

int
main(int argc,
     char **argv){
  string str1,str2,str3,appstr;
  DNA_b128_String dna1,dna2;
  simple_string_distance d1,d2;
  TN_string_distance t1,t2;
  
  float tmpTS,tmpTV,diffTS,diffTV;
  
  int test = 1;

  
  
  //  exit(1);
  //-----------------------
  cout <<"--------------\n Test " << (test++)<< endl <<
    "Simple test" <<endl;
  str1 = "aaaaaa";
  dna1 = DNA_b128_String(str1.length(),str1);
  dna2 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );
  assert( dna1 == dna2 );

  cout <<"--------------\n Test " << (test++)<<  endl <<
    "Simple test" <<endl;
  str1 = "cgcgcg";
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );
  
  cout <<"--------------\n Test " << (test++)<< endl <<
    "test 16 chars" <<endl;
  str1 = "cgcgcgcgcgcgcgcg";
  assert(str1.length() == 16);
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );
  
  cout <<"--------------\n Test " << (test++)<< endl <<
    "test 17 chars" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "c";
  assert(str1.length() == 17);
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );

  
  cout <<"--------------\n Test " << (test++)<< endl <<
    "test 32 chars" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg";
  assert(str1.length() == 32);
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );

  
  cout <<"--------------\n Test " << (test++)<< endl <<
    "test 33 chars" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "c";
  assert(str1.length() == 33);
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );


  cout <<"--------------\n Test " << (test++)<< endl <<
    "test 48 chars" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg";
  assert(str1.length() == 48);
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );

  
  cout <<"--------------\n Test " << (test++)<< endl <<
    "test 49 chars" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "c";
  assert(str1.length() == 49);
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );

  
  cout <<"--------------\n Test " << (test++)<< endl <<
    "test 64 chars" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg";
  assert(str1.length() == 64);
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );
  
  cout <<"--------------\n Test " << (test++)<< endl <<
    "test 65 chars" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "t";
  assert(str1.length() == 65);
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );


  cout <<"--------------\n Test " << (test++)<< endl <<
    "apppend test" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg";
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );
  appstr = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg";
  str1 += appstr;
  dna1.append(appstr);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );

  appstr = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "t";
  str1 += appstr;
  dna1.append(appstr);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );
  
  cout <<"--------------\n Test " << (test++)<< endl <<
    "apppend test with NOT_ALLOWED" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg";
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );
  appstr = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcZcgcg" "t";
  str1.append(appstr,0,appstr.find('Z'));
  dna1.append(appstr);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );
  
  cout <<"--------------\n Test " << (test++)<< endl <<
    "apppend maxRead test" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg";
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );
  appstr = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "t";
  str1.append(appstr);
  dna1.append(appstr);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );
  
  cout <<"--------------\n Test " << (test++)<< endl <<
    "apppend UNKNOWN test" <<endl;
  str1 = "cgcgc--gcgcgcgcg" "cgcgcgcgcgcgcgcg";
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );
  appstr = "cg-gc-cgcgcgcgcg" "cgcgcgc-cgcgcgcg" "t";
  dna1.append(appstr);
  str1 += appstr;
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );
    
  cout <<"--------------\n Test " << (test++)<< endl <<
    "apppend AMBIGIOUS test" <<endl;
  str1 = "cgcgc--gcgmgcgcg" "cgcgmgcgcgcgcgcg";
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );
  appstr = "cg-gc-cgcgngcgcg" "cgcgngc-cgcgcgcg" "t";
  dna1.append(appstr);
  str1 += appstr;
  appstr = "cg-gc-cgcgngcgcg" "cgcgngc-cgcgcgcg" "mrwsykvhdbn";
  dna1.append(appstr);
  str1 += appstr;
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );

  cout <<"--------------\n Test " << (test++)<< endl <<
    "set test" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgc-cgcg" "c-cgagcgtgcgcgcg" "cgcgcgcgcgcgcgcg" "c";
  assert(str1.length() == 65);
  dna1 = DNA_b128_String(str1.length(),str1);
  //dna1.setNucleotide(0,DNA_A_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );

  str1[0] ='t';dna1.setNucleotide(0,DNA_T_);
  str1[17]='t';dna1.setNucleotide(17,DNA_T_);
  str1[64]='t';dna1.setNucleotide(64,DNA_T_);
  str1[32]='t';dna1.setNucleotide(32,DNA_T_);
  str2 = dna1.toString();
  cout << str2 << endl;
  cout << str1 << endl;
  assert( str1 == str2 );

  cout <<"--------------\n Test " << (test++)<< endl <<
    "DIST TEST level 1  ... 64" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg";
  cout <<"creating strings of length " << str1.length() << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  str2 = str1;
  //A->G is a transition A->C a transversion
  str1[10] = 'a';//tv
  str1[8] = 'a';//tv
  str1[7] = 'a';//ts
  str1[0]= '-';//tv
  str1[63] = 'a';//ts
  dna2 = DNA_b128_String(str1.length(),str1);
  cout << dna1 << endl;
  cout << dna2 << endl;
  TEST_STRING_DISTANCE(dna1,dna2,str2,str1);

  cout <<"--------------\n Test " << (test++)<< endl <<
    "DIST TEST level 1 ... 128" <<endl;
  str1 = "cgcgcgcgcg-gcgcg" "cgc-cgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg";
  str1 = str1 + str1;
  cout <<"creating strings of length " << str1.length() << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  str2 = str1;
  //A->G is a transition A->C a transversion
  str1[10] = 'a';
  str1[8] = 'a';
  str1[7] = 'a';
  str1[0]= '-';
  str1[63] = 'a';
  dna2 = DNA_b128_String(str1.length(),str1);
  cout << dna1 << endl;
  cout << dna2 << endl;
  TEST_STRING_DISTANCE(dna1,dna2,str2,str1);    
  
  cout <<"--------------\n Test " << (test++)<< endl <<
    "DIST TEST level 2" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg";
  str1 = str1+str1+str1+str1;
  str1 = str1+str1+str1+str1;//i.e. sixteen b128
  cout <<"creating strings of length " << str1.length() << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  str2 = str1;
  //A->G is a transition A->C a transversion
  str1[10] = 'a';
  str1[8] = 'a';
  str1[7] = 'a';
  str1[0]= '-';
  str1[64] = 'a';
  str1[128] = 't';
  dna2 = DNA_b128_String(str1.length(),str1);
  TEST_STRING_DISTANCE(dna1,dna2,str2,str1);

  cout <<"--------------\n Test " << (test++)<< endl <<
    "DIST TEST level 3" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg";//one b128
  str1 = str1+str1+str1+str1;// four b128
  str1 = str1+str1+str1+str1;//i.e. sixteen b128
  str1 = str1+str1+str1+str1;//64 b128s
  str1 = str1+str1+str1+str1;//256 b128s
  
  cout <<"creating strings of length " << str1.length() << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  str2 = str1;
  //A->G is a transition A->C a transversion
  str1[10] = 'a';
  str1[8] = 'a';
  str1[7] = 'a';
  str1[0]= '-';
  str1[64] = 'a';
  str1[128] = 't';
  dna2 = DNA_b128_String(str1.length(),str1);
  TEST_STRING_DISTANCE(dna1,dna2,str2,str1);

  cout <<"--------------\n Test " << (test++)<< endl <<
    "DIST TEST level 4" <<endl;
  str1 = "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg" "cgcgcgcgcgcgcgcg";//one b128
  str1 = str1+str1+str1+str1;// four b128
  str1 = str1+str1+str1+str1;//i.e. sixteen b128
  str1 = str1+str1+str1+str1;//64 b128s
  str1 = str1+str1+str1+str1;//256 b128s
  str1 = str1+str1+str1+str1;//1024 b128s
  str1 = str1+str1+str1+str1;//4096 b128s
  str1 = str1+str1+str1+str1;//16384 b128s
  str1 = str1+str1+str1+str1;//... b128s
  
  cout <<"creating strings of length " << str1.length() << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  str2 = str1;
  //A->G is a transition A->C a transversion
  for ( unsigned int i = 0 ; i < str1.length()/2 ; i++ ){
    int rand_index = (int) (((float)str2.length())*rand()/(RAND_MAX+1.0));
    str2[rand_index] = 'a';
    if ( i % 5 == 0 )
      str2[rand_index] ='-';
  }
  dna2 = DNA_b128_String(str2.length(),str2);
  TEST_STRING_DISTANCE(dna1,dna2,str1,str2);
  
  cout <<"--------------\n Test " << (test++)<< endl <<
    "DIST TEST level 4" <<endl;
  str1 = "cgcgcattgcgattacgcgcg" "ttcgcgcgatatacgcgcgatatcgcg" "cgcgatatatcgatatcgcgcgcgcg" "cgcgcgcata";
  str1 = str1+str1+str1+str1;
  str1 = str1+str1+str1+str1;
  str1 = str1+str1+str1+str1;
  str1 = str1+str1+str1+str1;
  str1 = str1+str1+str1+str1;
  str1 = str1+str1+str1+str1;
  str1 = str1+str1+str1+str1;
  str1 = str1+str1+str1+str1;
  
  cout <<"creating strings of length " << str1.length() << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  str2 = str1;
  //A->G is a transition A->C a transversion
  for ( unsigned int i = 0 ; i < str1.length()/2 ; i++ ){
    int rand_index = (int) (((float)str1.length())*rand()/(RAND_MAX+1.0));
    str1[rand_index] = 'a';
    if ( i % 5 == 0 )
      str1[rand_index] ='-';
  }
  dna2 = DNA_b128_String(str1.length(),str1);
  TEST_STRING_DISTANCE(dna1,dna2,str2,str1);

  
  cout <<"--------------\n Test " << (test++)<< endl <<
    "DETAILED TEST" <<endl;
  str1 = "aaaaccccggggtttt";
  str2 = "acgtacgtacgtacgt";

  cout << str1 << endl;
  cout << str2 << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  dna2 = DNA_b128_String(str2.length(),str2);

  TEST_STRING_DISTANCE(dna1,dna2,str1,str2);


  //-----------------------
  // random length tests
  int numRandTests = 100;
  int maxLength = 20000;
  cout <<"--------------\n Random length Test "<< endl;

  for ( int i = 0 ; i < numRandTests ; i++){
    int length = (int) (((float)maxLength)*rand()/(RAND_MAX+1.0));
    str1.clear();
    cout << "-------------------\nLength="<< length <<"    capacity="<<str1.capacity() <<endl;
    for ( int j = 0 ; j < length ; j++ ){
      nucleotide n = (nucleotide)((int) (((float)4)*rand()/(RAND_MAX+1.0)));
      str1.append(1,nucleotide2char(n));
      if ( ((int)(((float)20)*rand()/(RAND_MAX+1.0))) == 1 )
        str1[str1.length()-1] ='-';
    }
    cout << "strlen=" << str1.length() << endl;
    dna1 = DNA_b128_String(str1.length(),str1);
    dna2 = DNA_b128_String(str1.length(),str1);
    assert( dna1 == dna2 );
    str2 = str1;
    //A->G is a transition A->C a transversion
    for ( unsigned int i = 0 ; i < str1.length()/2 ; i++ ){
      int rand_index = (int) (((float)str1.length())*rand()/(RAND_MAX+1.0));
      nucleotide n = (nucleotide)((int) (((float)4)*rand()/(RAND_MAX+1.0)));
      str1[rand_index] = nucleotide2char(n);
      if ( i % 5 == 0 )
        str1[rand_index] ='-';
    }
    dna2 = DNA_b128_String(str1.length(),str1);
    assert ( dna1 != dna2 );
    TEST_STRING_DISTANCE(dna1,dna2,str2,str1);
  }  

  //AMBIGUITY DISTANCE TEST
  cout <<"--------------\n AMBIGUITY DISTANCE"<< endl;
    
  str1 = "agctagct-";
  dna1 = DNA_b128_String(str1.length(),str1);
  assert( str1 == dna1.toString() );
  str2 = str1;
  str2[0] = 'c';
  str2[1] ='t';
  str2[4] = 'g';
  str2[6] ='t';
  cout << str1 << endl;
  cout << str2 << endl;
  dna2 = DNA_b128_String(str2.length(),str2);
  assert ( str2 == dna2.toString() );
  dna1.calcAmbiguityProbabilitiesUNIFORM();
  dna2.calcAmbiguityProbabilitiesUNIFORM();
  d1 = DNA_b128_String::computeDistance(dna1,dna2);
  d1 = DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(d1,dna1,dna2);
  cout << d1 << endl;

  cout <<"adding ambigious " << endl;
  str1 += "m";
  str2 += "a";
  cout << str1 << endl;
  cout << str2 << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  assert( str1 == dna1.toString() );
  dna2 = DNA_b128_String(str2.length(),str2);
  assert ( str2 == dna2.toString() );
  dna1.calcAmbiguityProbabilitiesUNIFORM();
  dna2.calcAmbiguityProbabilitiesUNIFORM();
  d1 = DNA_b128_String::computeDistance(dna1,dna2);
  d1 = DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(d1,dna1,dna2);
  cout << d1 << endl;
  cout <<"adding  ambigious " << endl;
  str1 += "n";
  str2 += "a";
  cout << str1 << endl;
  cout << str2 << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  assert( str1 == dna1.toString() );
  dna2 = DNA_b128_String(str2.length(),str2);
  assert ( str2 == dna2.toString() );
  dna1.calcAmbiguityProbabilitiesUNIFORM();
  dna2.calcAmbiguityProbabilitiesUNIFORM();
  d1 = DNA_b128_String::computeDistance(dna1,dna2);
  d1 = DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(d1,dna1,dna2);
  cout << d1 << endl;
  cout <<"adding  ambigious " << endl;
  str1 += "n";
  str2 += "n";
  cout << str1 << endl;
  cout << str2 << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  assert( str1 == dna1.toString() );
  dna2 = DNA_b128_String(str2.length(),str2);
  assert ( str2 == dna2.toString() );
  dna1.calcAmbiguityProbabilitiesUNIFORM();
  dna2.calcAmbiguityProbabilitiesUNIFORM();
  d1 = DNA_b128_String::computeDistance(dna1,dna2);
  d1 = DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(d1,dna1,dna2);
  cout << d1 << endl;

  cout <<"adding  ambigious " << endl;
  str1 += "mrwsykvhdbn";
  str2 += "aaaaaaaaaaa";
  cout << str1 << endl;
  cout << str2 << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  assert( str1 == dna1.toString() );
  dna2 = DNA_b128_String(str2.length(),str2);
  assert ( str2 == dna2.toString() );
  dna1.calcAmbiguityProbabilitiesUNIFORM();
  dna2.calcAmbiguityProbabilitiesUNIFORM();
  tmpTS = d1.transitions;
  tmpTV = d1.transversions;
  d1 = DNA_b128_String::computeDistance(dna1,dna2);
  d1 = DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(d1,dna1,dna2);
  diffTS = d1.transitions - tmpTS;
  diffTV = d1.transversions - tmpTV;  
  cout << d1 << endl;

  cout <<"adding  ambigious " << endl;
  str1 += "mrwsykvhdbn";
  str2 += "ccccccccccc";
  cout << str1 << endl;
  cout << str2 << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  assert( str1 == dna1.toString() );
  dna2 = DNA_b128_String(str2.length(),str2);
  assert ( str2 == dna2.toString() );
  dna1.calcAmbiguityProbabilitiesUNIFORM();
  dna2.calcAmbiguityProbabilitiesUNIFORM();
  tmpTS = d1.transitions;
  tmpTV = d1.transversions;
  d1 = DNA_b128_String::computeDistance(dna1,dna2);
  d1 = DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(d1,dna1,dna2);
  assert ( fabs(diffTS - (d1.transitions - tmpTS)) <= 0.0001 );
  assert ( fabs(diffTV - (d1.transversions - tmpTV)) <= 0.0001 );
  diffTS = d1.transitions - tmpTS;
  diffTV = d1.transversions - tmpTV;
  cout << d1 << endl;

  cout <<"adding  ambigious " << endl;
  str1 += "mrwsykvhdbn";
  str2 += "ggggggggggg";
  cout << str1 << endl;
  cout << str2 << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  assert( str1 == dna1.toString() );
  dna2 = DNA_b128_String(str2.length(),str2);
  assert ( str2 == dna2.toString() );
  dna1.calcAmbiguityProbabilitiesUNIFORM();
  dna2.calcAmbiguityProbabilitiesUNIFORM();
  tmpTS = d1.transitions;
  tmpTV = d1.transversions;
  d1 = DNA_b128_String::computeDistance(dna1,dna2);
  d1 = DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(d1,dna1,dna2);
  assert ( fabs(diffTS - (d1.transitions - tmpTS)) <= 0.0001 );
  assert ( fabs(diffTV - (d1.transversions - tmpTV)) <= 0.0001 );
  diffTS = d1.transitions - tmpTS;
  diffTV = d1.transversions - tmpTV;
  cout << d1 << endl;

  
  cout <<"adding  ambigious " << endl;
  str1 += "mrwsykvhdbn";
  str2 += "ttttttttttt";
  cout << str1 << endl;
  cout << str2 << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  assert( str1 == dna1.toString() );
  dna2 = DNA_b128_String(str2.length(),str2);
  assert ( str2 == dna2.toString() );
  dna1.calcAmbiguityProbabilitiesUNIFORM();
  dna2.calcAmbiguityProbabilitiesUNIFORM();
  tmpTS = d1.transitions;
  tmpTV = d1.transversions;
  d1 = DNA_b128_String::computeDistance(dna1,dna2);
  d1 = DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(d1,dna1,dna2);
  assert ( fabs(diffTS - (d1.transitions - tmpTS)) <= 0.0001 );
  assert ( fabs(diffTV - (d1.transversions - tmpTV)) <= 0.0001 );
  diffTS = d1.transitions - tmpTS;
  diffTV = d1.transversions - tmpTV;
  cout << d1 << endl;

  cout <<"adding  ambigious " << endl;
  str1 += "mrwsykvhdbn";
  str2 += "-----------";
  cout << str1 << endl;
  cout << str2 << endl;
  dna1 = DNA_b128_String(str1.length(),str1);
  assert( str1 == dna1.toString() );
  dna2 = DNA_b128_String(str2.length(),str2);
  assert ( str2 == dna2.toString() );
  dna1.calcAmbiguityProbabilitiesUNIFORM();
  dna2.calcAmbiguityProbabilitiesUNIFORM();
  d2 = d1;
  d1 = DNA_b128_String::computeDistance(dna1,dna2);
  d1 = DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(d1,dna1,dna2);
  cout <<"d2: " <<  d2 << endl;
  cout <<"d1: " <<  d1 << endl;
  assert ( d1.transversions == d2.transversions );
  assert ( d1.transitions == d2.transitions );
  assert ( d1.deletedPositions == d2.deletedPositions+11 ) ;
  cout << d1 << endl;

  SEPARATOR();SEPARATOR();SEPARATOR();
  timing_test();
    
  return 1;
}


// void
// parsi_test(){


//   string str1,str2,str3,appstr;
//   DNA_b128_String dna1,dna2;

//   //------------------------------------------------------
//   //PARSIMONY TEST
//   DNA_b128_String parsi;
//   cout <<"--------------\n PARSIMONY TEST"<< endl;
//   str1 = "a";
//   dna1 = DNA_b128_String(str1.length(),str1);
//   str2 = "a";
//   dna2 = DNA_b128_String(str2.length(),str2);
//   DNA_b128_String::create_weighted_parsimonious(parsi,dna1,dna2);

//   PRINT(dna1);
//   PRINT(dna2);
//   PRINT(parsi);
//   ASSERT_EQ(dna1,parsi);

//   cout <<"--------------\n";
//   str1 = "a";
//   dna1 = DNA_b128_String(str1.length(),str1);
//   str2 = "c";
//   dna2 = DNA_b128_String(str2.length(),str2);
//   DNA_b128_String::create_weighted_parsimonious(parsi,dna1,dna2);

//   PRINT(dna1);
//   PRINT(dna2);
//   PRINT(parsi);
//   assert(parsi.getNucleotide(0) == DNA_M_ );
  
  
//   cout <<"--------------\n";
//   str1 = "aaaaccccggggtttt";
//   dna1 = DNA_b128_String(str1.length(),str1);
//   str2 = "acgtacgtacgtacgt";
//   dna2 = DNA_b128_String(str2.length(),str2);
//   DNA_b128_String::create_weighted_parsimonious(parsi,dna1,dna2);

//   PRINT(dna1);
//   PRINT(dna2);
//   PRINT(parsi);
//   //the parsimony
//   str2 = "amrw" "mcsy" "rsgk" "wykt";
//   dna2 = DNA_b128_String(str2.length(),str2);
//   ASSERT_EQ(dna2,parsi);
  

//   cout <<"--------------\n";
//   str1 = "aaaaccccggggtttt" "mrws" "ykvh" "dbn"  "mrws" "ykvh" "dbn" "mrws" "ykvh" "dbn" "mrws" "ykvh" "dbn"
//     "mrws" "ykvh" "dbn"
//     "mrws" "ykvh" "dbn"
//     "mrws";
//   dna1 = DNA_b128_String(str1.length(),str1);
//   str2 = "acgtacgtacgtacgt" "aaaa" "aaaa" "aaa"  "cccc" "cccc" "ccc" "gggg" "gggg" "ggg" "tttt" "tttt" "ttt"
//     "rwsm" "mrws" "mrw"
//     "mrws" "ykvh" "dbn"
//     "vdhb";
//   dna2 = DNA_b128_String(str2.length(),str2);
//   DNA_b128_String::create_weighted_parsimonious(parsi,dna1,dna2);

//   PRINT(dna1);
//   PRINT(dna2);
//   PRINT(parsi);
//   //the parsimony
//   str2 = "amrw" "mcsy" "rsgk" "wykt"
//     "aaav" "hdaa" "ana"
//     "cvhc" "cbcc" "ncc"
//     "vgdg" "bggn" "ggg"
//     "hdtb" "ttnt" "ttt";
//   str2+= "aanc" "cgac" "agw"
//     "mrws" "ykvh" "dbn"
//     "mrws";
//   dna2 = DNA_b128_String(str2.length(),str2);
//   PRINT(dna2);
//   ASSERT_EQ(dna2,parsi);
//   parsi.printVerbose(cout);

// }

