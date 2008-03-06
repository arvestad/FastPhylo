//--------------------------------------------------
//                                        
// File: computeDistance_DNA_b128_String.cpp
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: new_comp.cpp,v 1.3 2006/12/20 14:58:03 isaac Exp $                                 
//
//--------------------------------------------------

#include <string>

#include "DNA_b128_String.hpp"
#include <iostream>

using namespace std;
//----------------------
// LEVEL SUMS etc
#define DEBUG_MAIN(INP)  if(false){INP;}


//------------------------------------
// CONVERTING LEVELS
// This macro takes three arguments:
// SUM_CURRENT  -  The b128 sum of the current level for either TS or TV.
// SUM_PREVIOUS -  The b128 sum of the previous level for either TS or TV.
// MASK         -  A mask one ones in all positions of every other block of the previous level.
// SHIFT        -  The number of positions to right shift the sum of the previous level to get them above each other.
//
// Example usage:
// SUM_WITH_PREVIOUS_LEVEL(sum_ts_l2, sum_ts_l1, TWO_BIT_MASK, TWO);
// Takes the sums from level 1 and adds them to the approraite blocks of sum_ts_l2.

#define SUM_WITH_PREVIOUS_LEVEL(SUM_CURRENT,SUM_PREVIOUS,MASK,SHIFT)\
{SUM_CURRENT = add_b128(SUM_CURRENT,add_b128(and_b128(SUM_PREVIOUS,MASK),and_b128(shift_each32_bits_right_b128(SUM_PREVIOUS,SHIFT),MASK)));}

//the same but shifts whole bytes instead of bits it is probably
//faster than the previous one since it utilzes immediates instead of
//variables.
#define SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(SUM_CURRENT,SUM_PREVIOUS,MASK,BYTESHIFT)\
{SUM_CURRENT = add_b128(SUM_CURRENT,add_b128(and_b128(SUM_PREVIOUS,MASK),and_b128(shift_bytes_right_b128(SUM_PREVIOUS,BYTESHIFT),MASK)));}

// CONVERTS THE SUM
// Takes the sum of the previous level adds all two adjacent blocks into
// one block of the next level.
#define CONVERT_SUM(SUM_CURRENT,SUM_PREVIOUS,MASK,SHIFT)\
{SUM_CURRENT = add_b128(and_b128(SUM_PREVIOUS,MASK),and_b128(shift_each32_bits_right_b128(SUM_PREVIOUS,SHIFT),MASK));}

//the same but shifts whole bytes instead of bits
#define CONVERT_SUM_IMMEDIATEBYTESHIFT(SUM_CURRENT,SUM_PREVIOUS,MASK,BYTESHIFT)\
{SUM_CURRENT = add_b128(and_b128(SUM_PREVIOUS,MASK),and_b128(shift_bytes_right_b128(SUM_PREVIOUS,BYTESHIFT),MASK));}



//------------------------------------
// DISTANCE COMPUTATION
simple_string_distance
DNA_b128_String::computeDistance(const DNA_b128_String &s1,
                                 const DNA_b128_String &s2){

  const b128 LEAST_SIGNIFCANT_BIT = set_all_bytes(0x55);
  const b128 ONE = set_first_int_b128(1);
  
  //LEVEL 2
  b128 TWO_BIT_MASK = set_all_bytes(0x33);
  const b128 TWO = set_first_int_b128(2);
  
  //LEVEL 3
  
  // A mask with every four bits set as follows "..00001111"
  const b128 FOUR_BIT_MASK = set_all_bytes(0x0f);
 const b128 FOUR = set_first_int_b128(4);
 
 //LEVEL 4
 const b128 EIGHT_BIT_MASK = set_all_shorts(0x00ff);
 const b128 EIGHT = set_first_int_b128(8);

 const b128 SIXTEEN_BIT_MASK = set_all_ints(0x0000ffff);
 const b128 SIXTEEN = set_first_int_b128(16);

  //The date that the functions will work on.
  //Note! it is only in this function and in dist_level_1()
  //that data is actually read.
  b128 *ptr1 = s1.data;
  b128 *ptr2 = s2.data;
  b128 *del_ptr1 = s1.unknownData;
  b128 *del_ptr2 = s2.unknownData;
  assert ( s1.getNumChars() == s2.getNumChars() );

  //----
  //Compute the number of times each level should be called.
  //
  // The general layout of the different levels is as follows:
  // dist_level_4() calls dist_level_3()
  // dist_level_3() calls dist_level_2()
  // dist_level_2() calls dist_level_1() 
  //
  // The block sizes at each level is as follows:
  // level 4: 16 bits, representing at most 65,535 missmatches
  // level 3: 8 bits, representing at most 255 missmatches
  // level 2: 4 bits, representing at most 15 missmatches
  // level 1: 2 bits, representing at most 3 missmatches
  //
  // At each level two blocks of the previous level are added into a
  // block. Thus each level can only call a specific number of times
  // to the other level before the blocks overflow. The number of
  // calls to the previous levels are:
  // level 4 -> level 3: 65,535/(255+255)=128
  // level 3 -> level 2: 255/(15+15)=8
  // level 2 -> level 1: 15/(3+3)=2
  //
  // Since dist_level_1() compairs three b128s after a call to
  // dist_level_X() a certain number of b128 have been read and
  // handled. The number of handled b128s for each level is:
  //level 4: 128*8*2*3 = 6144
  //level 3: 8*2*3     = 48
  //level 2: 2*3       = 6
  //level 1: 3         = 3  
  //
  // We now compute the number of times each level should be called:
  
  const int numDatasToRead = s1.getNumUsedDatas();
  const b128 *endOf_p1 = ptr1+numDatasToRead;

  //BIG LEVEL
  b128 total_sum_tv = set_zero_b128();
  b128 total_sum_ts = set_zero_b128();
  b128 total_sum_del = set_zero_b128();
  
  //level 4 
 while( true ){
    b128 sum_tv_l4 = set_zero_b128();
    b128 sum_ts_l4 = set_zero_b128();
    b128 sum_del_l4 = set_zero_b128();
    
    for(size_t l3iter=0; l3iter<128; l3iter++){
      b128 sum_tv_l3 = set_zero_b128();
      b128 sum_ts_l3 = set_zero_b128();
      b128 sum_del_l3 = set_zero_b128();

      for(size_t l2iter=0; l2iter<8; l2iter++){
	b128 sum_tv_l2 = set_zero_b128();
	b128 sum_ts_l2 = set_zero_b128();
	b128 sum_del_l2 = set_zero_b128();

	for(size_t l1iter=0; l1iter<2 ; l1iter++){
	  b128 del = or_b128(get_b128(del_ptr1),get_b128(del_ptr2));	  
	  b128 diff = andnot_b128(del,xor_b128(get_b128(ptr1),get_b128(ptr2)));
	  ptr1++;ptr2++;del_ptr1++; del_ptr2++;
	  _mm_prefetch((char*) del_ptr1,_MM_HINT_NTA);
	  _mm_prefetch((char*) del_ptr2,_MM_HINT_NTA);
	  _mm_prefetch((char*) ptr1,_MM_HINT_NTA);
	  _mm_prefetch((char*) ptr2,_MM_HINT_NTA);

	  b128 sum_tv_l1 = and_b128(shift_each32_bits_right_b128(diff,ONE),LEAST_SIGNIFCANT_BIT);
	  b128 sum_del_l1 = and_b128(del,LEAST_SIGNIFCANT_BIT);
	  b128 sum_ts_l1 = andnot_b128(sum_tv_l1, and_b128(diff,LEAST_SIGNIFCANT_BIT));
	  b128 tmp_tv;
  
	  if(endOf_p1==ptr1) goto SUM_UP_ALL_LEVELS_AND_BREAK;

	  del = or_b128(get_b128(del_ptr1),get_b128(del_ptr2));	  
	  diff = andnot_b128(del,xor_b128(get_b128(ptr1),get_b128(ptr2)));
	  ptr1++;ptr2++;del_ptr1++; del_ptr2++;
	  _mm_prefetch((char*) del_ptr1,_MM_HINT_NTA);
	  _mm_prefetch((char*) del_ptr2,_MM_HINT_NTA);
	  _mm_prefetch((char*) ptr1,_MM_HINT_NTA);
	  _mm_prefetch((char*) ptr2,_MM_HINT_NTA);

	  tmp_tv = and_b128(shift_each32_bits_right_b128(diff,ONE),LEAST_SIGNIFCANT_BIT);
	  sum_del_l1 = add_b128(sum_del_l1, and_b128(del,LEAST_SIGNIFCANT_BIT));
	  sum_tv_l1 = add_b128(sum_tv_l1, tmp_tv);
	  sum_ts_l1 = add_b128(sum_ts_l1,andnot_b128(tmp_tv, and_b128(diff,LEAST_SIGNIFCANT_BIT)));
	  
	  if(endOf_p1==ptr1) goto SUM_UP_ALL_LEVELS_AND_BREAK;
	  
	  del = or_b128(get_b128(del_ptr1),get_b128(del_ptr2));	  
	  diff = andnot_b128(del,xor_b128(get_b128(ptr1),get_b128(ptr2)));
	  ptr1++;ptr2++;del_ptr1++; del_ptr2++;
	  _mm_prefetch((char*) del_ptr1,_MM_HINT_NTA);
	  _mm_prefetch((char*) del_ptr2,_MM_HINT_NTA);
	  _mm_prefetch((char*) ptr1,_MM_HINT_NTA);
	  _mm_prefetch((char*) ptr2,_MM_HINT_NTA);
	  
	  tmp_tv = and_b128(shift_each32_bits_right_b128(diff,ONE),LEAST_SIGNIFCANT_BIT);
	  sum_del_l1 = add_b128(sum_del_l1, and_b128(del,LEAST_SIGNIFCANT_BIT));
	  sum_tv_l1 = add_b128(sum_tv_l1, tmp_tv);
	  sum_ts_l1 = add_b128(sum_ts_l1,andnot_b128(tmp_tv, and_b128(diff,LEAST_SIGNIFCANT_BIT)));
	  
	  SUM_WITH_PREVIOUS_LEVEL(sum_ts_l2,sum_ts_l1, TWO_BIT_MASK, TWO);
	  SUM_WITH_PREVIOUS_LEVEL(sum_tv_l2,sum_tv_l1, TWO_BIT_MASK, TWO);
	  SUM_WITH_PREVIOUS_LEVEL(sum_del_l2,sum_del_l1, TWO_BIT_MASK, TWO);

	  if(endOf_p1==ptr1) goto SUM_UP_ALL_LEVELS_AND_BREAK;
	  
	  continue;
	  
	  // WHEN ALL DATAS HAVE BEEN READ
	SUM_UP_ALL_LEVELS_AND_BREAK:
	  SUM_WITH_PREVIOUS_LEVEL(sum_ts_l2,sum_ts_l1, TWO_BIT_MASK, TWO);
	  SUM_WITH_PREVIOUS_LEVEL(sum_tv_l2,sum_tv_l1, TWO_BIT_MASK, TWO);
	  SUM_WITH_PREVIOUS_LEVEL(sum_del_l2,sum_del_l1, TWO_BIT_MASK, TWO);
	  
	  SUM_WITH_PREVIOUS_LEVEL(sum_ts_l3,sum_ts_l2, FOUR_BIT_MASK, FOUR);
	  SUM_WITH_PREVIOUS_LEVEL(sum_tv_l3,sum_tv_l2, FOUR_BIT_MASK, FOUR);
	  SUM_WITH_PREVIOUS_LEVEL(sum_del_l3,sum_del_l2, FOUR_BIT_MASK, FOUR);

	  SUM_WITH_PREVIOUS_LEVEL(sum_ts_l3,sum_ts_l2, FOUR_BIT_MASK, FOUR);
	  SUM_WITH_PREVIOUS_LEVEL(sum_tv_l3,sum_tv_l2, FOUR_BIT_MASK, FOUR);
	  SUM_WITH_PREVIOUS_LEVEL(sum_del_l3,sum_del_l2, FOUR_BIT_MASK, FOUR);
	  goto END_OF_READ_LOOP;
	}
	SUM_WITH_PREVIOUS_LEVEL(sum_ts_l3,sum_ts_l2, FOUR_BIT_MASK, FOUR);
	SUM_WITH_PREVIOUS_LEVEL(sum_tv_l3,sum_tv_l2, FOUR_BIT_MASK, FOUR);
	SUM_WITH_PREVIOUS_LEVEL(sum_del_l3,sum_del_l2, FOUR_BIT_MASK, FOUR);
      }
      SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(sum_ts_l4,sum_ts_l3, EIGHT_BIT_MASK, 1);
      SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(sum_tv_l4,sum_tv_l3, EIGHT_BIT_MASK, 1);
      SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(sum_del_l4,sum_del_l3, EIGHT_BIT_MASK, 1);
    }
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(total_sum_ts,sum_ts_l4, SIXTEEN_BIT_MASK, 2);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(total_sum_tv,sum_tv_l4, SIXTEEN_BIT_MASK, 2);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(total_sum_del,sum_del_l4, SIXTEEN_BIT_MASK,2);
  }
  
 END_OF_READ_LOOP:
  //-------------
  // COMBINE THE SUMS IN EACH INT
  //
  // Now each block in total_sum* is 32 bits. Thus to get the final
  // sum we simply need to add the four ints.

  total_sum_ts = add_b128(total_sum_ts, shift_bytes_right_b128(total_sum_ts,4));
  total_sum_ts = add_b128(total_sum_ts, shift_bytes_right_b128(total_sum_ts,8));
  total_sum_tv = add_b128(total_sum_tv, shift_bytes_right_b128(total_sum_tv,4));
  total_sum_tv = add_b128(total_sum_tv, shift_bytes_right_b128(total_sum_tv,8));
  total_sum_del = add_b128(total_sum_del, shift_bytes_right_b128(total_sum_del,4));
  total_sum_del = add_b128(total_sum_del, shift_bytes_right_b128(total_sum_del,8));
  
  simple_string_distance d= {get_int_0_b128(total_sum_del),
                   get_int_0_b128(total_sum_ts),
                   get_int_0_b128(total_sum_tv)};
  
  DEBUG_MAIN( cout << "MAIN_S = "<< d.transitions << endl);
  DEBUG_MAIN( cout << "MAIN_V = "<< d.transversions << endl);
  DEBUG_MAIN( cout << "MAIN_D = "<< d.deletedPositions << endl);
  
  return d;
}
