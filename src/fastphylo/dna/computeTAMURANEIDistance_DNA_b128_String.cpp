//--------------------------------------------------
//                                        
// File: computeDistance_DNA_b128_String.cpp
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: computeTAMURANEIDistance_DNA_b128_String.cpp,v 1.5 2006/12/10 19:58:53 isaac Exp $                        
//
//--------------------------------------------------

#include <string>

#include "fastphylo/dna/DNA_b128_String.hpp"
#include <iostream>

using namespace std;

//
// This file contains the implementation of:
// TN_string_distance
// DNA_b128_String::computeTAMURANEIDistance(const DNA_b128_String &s1,
//                                           const DNA_b128_String &s2)
//
// The function computes the number of purine, pyrimidine
// transitions and transversions between two strings.
//
// See computeDistance_DNA_b128_String.cpp for further documentation.


//-----------------------------
// Declarations of methods that compute the sums
// at different levels
static void dist_level_1(b128 &sum_ts_l1, b128 &sum_pyrts_l1, b128 &sum_tv_l1, b128 &sum_del_l1);
static void dist_level_2(b128 &sum_ts_l2, b128 &sum_pyrts_l2, b128 &sum_tv_l2, b128 &sum_del_l2);
static void dist_level_3(b128 &sum_ts_l3, b128 &sum_pyrts_l3, b128 &sum_tv_l3, b128 &sum_del_l3);
static void dist_level_4(b128 &sum_ts_l4, b128 &sum_pyrts_l4, b128 &sum_tv_l4, b128 &sum_del_l4);



//-------------------------
// THE DATA
//
// All functions below work on this memory. The pointers are
// initialized in computeDistance().

static b128 *ptr1;
static b128 *ptr2;

static b128 *del_ptr1;
static b128 *del_ptr2;

static void inline 
PREFETCH_DATA(size_t numDatas){
  b128 *_del_ptr1;
  b128 *_del_ptr2;
  b128 *_ptr1;
  b128 *_ptr2;
  
  _del_ptr1 = del_ptr1;
  _del_ptr2 = del_ptr2;
  _ptr1 = ptr1;
  _ptr2 = ptr2;
  
  for(size_t i = 0 ; i<numDatas ; i++){
    _mm_prefetch((char*) (_del_ptr1++),_MM_HINT_NTA);
    _mm_prefetch((char*) (_del_ptr2++),_MM_HINT_NTA);
    _mm_prefetch((char*) (_ptr1++),_MM_HINT_NTA);
    _mm_prefetch((char*) (_ptr2++),_MM_HINT_NTA);
  }
}

//----------------------
// LEVEL SUMS etc

// sse2_wrapper.h's set_*_b128()/set_first_int_b128() are extern "C"
// wrappers around SSE2/simde intrinsics - they cannot throw, but being
// plain C they carry no noexcept for clang-tidy to see.
// NOLINTBEGIN(bugprone-throwing-static-initialization)

// LEVEL 1
//static b128 sum_ts_l1;
//static b128 sum_pyrts_l1;
//static b128 sum_tv_l1;
//static b128 sum_del_l1;
// A mask where the least significant bit in every two block is set.
static b128 LEAST_SIGNIFCANT_BIT = set_all_bytes(0x55);
static b128 ONE = set_first_int_b128(1);

//LEVEL 2
//static b128 sum_ts_l2;
//static b128 sum_pyrts_l2;
//static b128 sum_tv_l2;
//static b128 sum_del_l2;
// A mask with every two bits set as follows "..0011011"
static b128 TWO_BIT_MASK = set_all_bytes(0x33);
static b128 TWO = set_first_int_b128(2);

//LEVEL 3
//static b128 sum_ts_l3;
//static b128 sum_pyrts_l3;
//static b128 sum_tv_l3;
//static b128 sum_del_l3;
// A mask with every four bits set as follows "..00001111"
static b128 FOUR_BIT_MASK = set_all_bytes(0x0f);
static b128 FOUR = set_first_int_b128(4);

//LEVEL 4
//static b128 sum_ts_l4;
//static b128 sum_pyrts_l4;
//static b128 sum_tv_l4;
//static b128 sum_del_l4;
// A mask with every eigth bits set as follows "..0000 0000 1111 1111"
static b128 EIGHT_BIT_MASK = set_all_shorts(0x00ff);
static b128 EIGHT = set_first_int_b128(8);

//FINAL LEVEL
static b128 SIXTEEN_BIT_MASK = set_all_ints(0x0000ffff);
static b128 SIXTEEN = set_first_int_b128(16);
// NOLINTEND(bugprone-throwing-static-initialization)

//------------------------------------
// CONVERTING LEVELS
//
// sum_current    -  The b128 sum of the current level for either TS or TV.
// sum_previous   -  The b128 sum of the previous level for either TS or TV.
// mask           -  A mask with ones in all positions of every other block of the previous level.
// shift          -  The number of positions to right shift the sum of the previous level to get them above each other.
//
// Example usage:
// sum_with_previous_level(sum_ts_l2, sum_ts_l1, TWO_BIT_MASK, TWO);
// Takes the sums from level 1 and adds them to the approraite blocks of sum_ts_l2.
static inline void
sum_with_previous_level(b128 &sum_current, b128 sum_previous, b128 mask, b128 shift){
  sum_current = add_b128(sum_current,add_b128(and_b128(sum_previous,mask),and_b128(shift_each32_bits_right_b128(sum_previous,shift),mask)));
}

//the same but shifts whole bytes instead of bits, probably faster than
//the previous one since it utilizes immediates instead of variables.
//shift_bytes_right_b128<N>()'s byte count is a template parameter, not
//a macro parameter (see sse2_wrapper.h) - still needs the surrounding
//braces (bugprone-macro-parentheses doesn't apply to the template
//argument itself, only to SUM_CURRENT/SUM_PREVIOUS/MASK below).
#define SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(SUM_CURRENT,SUM_PREVIOUS,MASK,BYTESHIFT)\
{(SUM_CURRENT) = add_b128((SUM_CURRENT),add_b128(and_b128((SUM_PREVIOUS),(MASK)),and_b128(shift_bytes_right_b128<(BYTESHIFT)>((SUM_PREVIOUS)),(MASK))));}

// CONVERTS THE SUM
// Takes the sum of the previous level adds all two adjacent blocks into
// one block of the next level.
static inline void
convert_sum(b128 &sum_current, b128 sum_previous, b128 mask, b128 shift){
  sum_current = add_b128(and_b128(sum_previous,mask),and_b128(shift_each32_bits_right_b128(sum_previous,shift),mask));
}

//the same but shifts whole bytes instead of bits - same
//shift_bytes_right_b128<N>() template as the sibling above.
#define CONVERT_SUM_IMMEDIATEBYTESHIFT(SUM_CURRENT,SUM_PREVIOUS,MASK,BYTESHIFT)\
{(SUM_CURRENT) = add_b128(and_b128((SUM_PREVIOUS),(MASK)),and_b128(shift_bytes_right_b128<(BYTESHIFT)>((SUM_PREVIOUS)),(MASK)));}



//------------------------------------
// REMAINDER HANDLING
//
// dist_level_1() always consumes three b128s at a time, so at most two
// are left over once the level_1/2/3/4 call counts are computed. Sums
// them directly into total_sum_ts/pyrts/tv/del - a fallthrough switch
// by design: case 2 handles the first leftover (there's nothing summed
// yet, so it can assign rather than add), then falls into case 1's
// body, which handles the second leftover (or the only one, if there
// was just one) by adding into whatever case 2 already produced. This
// was the single largest contributor to computeTAMURANEIDistance()'s
// cognitive complexity (85 of it, on its own) - extracted here
// verbatim (down to the asserts and pointer increments) rather than
// restructured, since the fallthrough *is* the logic, not an
// accident: rewriting it as two independent branches would either
// duplicate the case-1 body or require an extra bool parameter.
static void
sumRemainderBlocks(int rest_num_b128s, b128 &total_sum_ts, b128 &total_sum_pyrts, b128 &total_sum_tv, b128 &total_sum_del){
  b128 diff;
  b128 del;
  b128 tmp_tv;
  switch( rest_num_b128s ){
  case 2:
    assert ( equal_b128(total_sum_tv,set_zero_b128()) );//assuming that nothing summed so far.
    assert ( equal_b128(total_sum_ts,set_zero_b128()) );//assuming that nothing summed so far.
    assert ( equal_b128(total_sum_pyrts,set_zero_b128()) );//assuming that nothing summed so far.
    assert ( equal_b128(total_sum_del,set_zero_b128()) );//assuming that nothing summed so far.

    del = or_b128(get_b128(del_ptr1),get_b128(del_ptr2));
    diff = andnot_b128(del,xor_b128(get_b128(ptr1),get_b128(ptr2)));

    total_sum_del = and_b128(del,LEAST_SIGNIFCANT_BIT);

    tmp_tv = and_b128(shift_each32_bits_right_b128(diff,ONE),LEAST_SIGNIFCANT_BIT);
    total_sum_tv = tmp_tv;
    total_sum_ts = andnot_b128(tmp_tv, and_b128(diff,LEAST_SIGNIFCANT_BIT));
    total_sum_pyrts = and_b128(shift_each32_bits_right_b128(and_b128(get_b128(ptr1),get_b128(ptr2)),ONE),total_sum_ts);

    ++ptr1;++ptr2;
    ++del_ptr1;++del_ptr2;
    [[fallthrough]];
  case 1: {
    del = or_b128(get_b128(del_ptr1),get_b128(del_ptr2));

    diff = andnot_b128(del,xor_b128(get_b128(ptr1),get_b128(ptr2)));

    total_sum_del = add_b128(total_sum_del, and_b128(del,LEAST_SIGNIFCANT_BIT));

    tmp_tv = and_b128(shift_each32_bits_right_b128(diff,ONE),LEAST_SIGNIFCANT_BIT);
    total_sum_tv = add_b128(total_sum_tv, tmp_tv);
    b128 tmp_ts = andnot_b128(tmp_tv, and_b128(diff,LEAST_SIGNIFCANT_BIT));
    total_sum_ts = add_b128(total_sum_ts,tmp_ts);
    total_sum_pyrts = add_b128(total_sum_pyrts,and_b128(shift_each32_bits_right_b128(and_b128(get_b128(ptr1),get_b128(ptr2)),ONE),tmp_ts));

    ++ptr1;++ptr2;
    ++del_ptr1;++del_ptr2;
    break;
  }
  case 0:
    break; // nothing left to sum
  default:
    assert(false && "rest_num_b128s must be 0, 1 or 2");
    break;
  }
}

//------------------------------------
// DISTANCE COMPUTATION
TN_string_distance
DNA_b128_String::computeTAMURANEIDistance(const DNA_b128_String &s1,
                                 const DNA_b128_String &s2){

  //The date that the functions will work on.
  //Note! it is only in this function and in dist_level_1()
  //that data is actually read.
  ptr1 = s1.data;
  ptr2 = s2.data;
  del_ptr1 = s1.unknownData;
  del_ptr2 = s2.unknownData;
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
  
  int rest_num_b128s = s1.getNumUsedDatas();

  int num_level_4 = rest_num_b128s / 6144; 
  rest_num_b128s = rest_num_b128s % 6144;  
  
  int num_level_3 = rest_num_b128s / 48; 
  rest_num_b128s = rest_num_b128s % 48;  

  int num_level_2 = rest_num_b128s / 6; 
  rest_num_b128s = rest_num_b128s % 6;  

  int num_level_1 = rest_num_b128s / 3; 
  rest_num_b128s = rest_num_b128s % 3;  

  //----
  // CALL THE LEVELS

  //this variable will contain the sums
  //the block size will be updated regularly
  b128 total_sum_tv = set_zero_b128();
  b128 total_sum_ts = set_zero_b128();
  b128 total_sum_pyrts = set_zero_b128();
  b128 total_sum_del = set_zero_b128();
  
  //Compute the remaining b128s. There are atmost two remaining.
  PREFETCH_DATA(rest_num_b128s);
  sumRemainderBlocks(rest_num_b128s, total_sum_ts, total_sum_pyrts, total_sum_tv, total_sum_del);



  //update the block size to size 4
  convert_sum(total_sum_ts,total_sum_ts,TWO_BIT_MASK,TWO);
  convert_sum(total_sum_pyrts,total_sum_pyrts,TWO_BIT_MASK,TWO);
  convert_sum(total_sum_tv,total_sum_tv,TWO_BIT_MASK,TWO);
  convert_sum(total_sum_del,total_sum_del,TWO_BIT_MASK,TWO);
  
  
  //level 1, num_level_1 is atmost 2
  //   switch ( num_level_1 ){
  //   case 2:
  //     dist_level_1();
  //     sum_with_previous_level(total_sum_ts,sum_ts_l1, TWO_BIT_MASK, TWO);
  //     sum_with_previous_level(total_sum_tv,sum_tv_l1, TWO_BIT_MASK, TWO);
  //     sum_with_previous_level(total_sum_del,sum_del_l1, TWO_BIT_MASK, TWO);
  //   case 1:
  //     dist_level_1();
  //     sum_with_previous_level(total_sum_ts,sum_ts_l1, TWO_BIT_MASK, TWO);
  //     sum_with_previous_level(total_sum_tv,sum_tv_l1, TWO_BIT_MASK, TWO);
  //     sum_with_previous_level(total_sum_del,sum_del_l1, TWO_BIT_MASK, TWO);
  //   }
  b128 sum_ts_l1;
  b128 sum_pyrts_l1;
  b128 sum_tv_l1;
  b128 sum_del_l1;
  PREFETCH_DATA(static_cast<size_t>(num_level_1)*3);
  for (  ; num_level_1 != 0 ; num_level_1-- ){    
    dist_level_1(sum_ts_l1,  sum_pyrts_l1, sum_tv_l1,  sum_del_l1);
    sum_with_previous_level(total_sum_ts,sum_ts_l1, TWO_BIT_MASK, TWO);
    sum_with_previous_level(total_sum_pyrts,sum_pyrts_l1, TWO_BIT_MASK, TWO);
    sum_with_previous_level(total_sum_tv,sum_tv_l1, TWO_BIT_MASK, TWO);
    sum_with_previous_level(total_sum_del,sum_del_l1, TWO_BIT_MASK, TWO);
  }

  
  //update the block size to size 8
  convert_sum(total_sum_ts,total_sum_ts,FOUR_BIT_MASK,FOUR);
  convert_sum(total_sum_pyrts,total_sum_pyrts,FOUR_BIT_MASK,FOUR);
  convert_sum(total_sum_tv,total_sum_tv,FOUR_BIT_MASK,FOUR);
  convert_sum(total_sum_del,total_sum_del,FOUR_BIT_MASK,FOUR);
  
  
  //level 2, num_level_2 is atmost 8
  b128 sum_ts_l2;
  b128 sum_pyrts_l2;
  b128 sum_tv_l2;
  b128 sum_del_l2;
  for ( ; num_level_2 != 0 ; num_level_2-- ){
    dist_level_2(sum_ts_l2,  sum_pyrts_l2, sum_tv_l2,  sum_del_l2);
    sum_with_previous_level(total_sum_ts,sum_ts_l2, FOUR_BIT_MASK, FOUR);
    sum_with_previous_level(total_sum_pyrts,sum_pyrts_l2, FOUR_BIT_MASK, FOUR);
    sum_with_previous_level(total_sum_tv,sum_tv_l2, FOUR_BIT_MASK, FOUR);
    sum_with_previous_level(total_sum_del,sum_del_l2, FOUR_BIT_MASK, FOUR);
  }

  
  //update the block size to size 16
  //  convert_sum(total_sum_ts,total_sum_ts,EIGHT_BIT_MASK,EIGHT);
  //  convert_sum(total_sum_tv,total_sum_tv,EIGHT_BIT_MASK,EIGHT);
  //  convert_sum(total_sum_del,total_sum_del,EIGHT_BIT_MASK,EIGHT);  
  CONVERT_SUM_IMMEDIATEBYTESHIFT(total_sum_ts,total_sum_ts,EIGHT_BIT_MASK,1);
  CONVERT_SUM_IMMEDIATEBYTESHIFT(total_sum_pyrts,total_sum_pyrts,EIGHT_BIT_MASK,1);
  CONVERT_SUM_IMMEDIATEBYTESHIFT(total_sum_tv,total_sum_tv,EIGHT_BIT_MASK,1);
  CONVERT_SUM_IMMEDIATEBYTESHIFT(total_sum_del,total_sum_del,EIGHT_BIT_MASK,1);
  
  
  //level 3, num_level_3 is atmost 128
  b128 sum_ts_l3;
  b128 sum_pyrts_l3;
  b128 sum_tv_l3;
  b128 sum_del_l3;
  for ( ; num_level_3 != 0 ; num_level_3-- ){
    dist_level_3(sum_ts_l3,  sum_pyrts_l3, sum_tv_l3,  sum_del_l3);
    //    sum_with_previous_level(total_sum_ts,sum_ts_l3, EIGHT_BIT_MASK, EIGHT);
    //    sum_with_previous_level(total_sum_tv,sum_tv_l3, EIGHT_BIT_MASK, EIGHT);
    //    sum_with_previous_level(total_sum_del,sum_del_l3, EIGHT_BIT_MASK, EIGHT);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(total_sum_ts,sum_ts_l3, EIGHT_BIT_MASK, 1);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(total_sum_pyrts,sum_pyrts_l3, EIGHT_BIT_MASK, 1);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(total_sum_tv,sum_tv_l3, EIGHT_BIT_MASK, 1);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(total_sum_del,sum_del_l3, EIGHT_BIT_MASK,1);        
  }
  
  
  //update the block size to size 32
  //  convert_sum(total_sum_ts,total_sum_ts, SIXTEEN_BIT_MASK,SIXTEEN);
  //  convert_sum(total_sum_tv,total_sum_tv, SIXTEEN_BIT_MASK,SIXTEEN);
  //  convert_sum(total_sum_del,total_sum_del, SIXTEEN_BIT_MASK,SIXTEEN);
  CONVERT_SUM_IMMEDIATEBYTESHIFT(total_sum_ts,total_sum_ts, SIXTEEN_BIT_MASK,2);
  CONVERT_SUM_IMMEDIATEBYTESHIFT(total_sum_pyrts,total_sum_pyrts, SIXTEEN_BIT_MASK,2);
  CONVERT_SUM_IMMEDIATEBYTESHIFT(total_sum_tv,total_sum_tv, SIXTEEN_BIT_MASK,2);
  CONVERT_SUM_IMMEDIATEBYTESHIFT(total_sum_del,total_sum_del, SIXTEEN_BIT_MASK,2);

  
  //level 4, num_level_4 is atmost Any number
  b128 sum_ts_l4;
  b128 sum_pyrts_l4;
  b128 sum_tv_l4;
  b128 sum_del_l4;
  for (  ; num_level_4 != 0 ; num_level_4-- ){
    dist_level_4(sum_ts_l4,  sum_pyrts_l4, sum_tv_l4,  sum_del_l4);
    //    sum_with_previous_level(total_sum_ts,sum_ts_l4, SIXTEEN_BIT_MASK, SIXTEEN);
    //    sum_with_previous_level(total_sum_tv,sum_tv_l4, SIXTEEN_BIT_MASK, SIXTEEN);
    //    sum_with_previous_level(total_sum_del,sum_del_l4, SIXTEEN_BIT_MASK, SIXTEEN);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(total_sum_ts,sum_ts_l4, SIXTEEN_BIT_MASK, 2);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(total_sum_pyrts,sum_pyrts_l4, SIXTEEN_BIT_MASK, 2);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(total_sum_tv,sum_tv_l4, SIXTEEN_BIT_MASK, 2);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(total_sum_del,sum_del_l4, SIXTEEN_BIT_MASK,2);

  }


  
  
  //-------------
  // COMBINE THE SUMS IN EACH INT
  //
  // Now each block in total_sum* is 32 bits. Thus to get the final
  // sum we simply need to add the four ints.

  total_sum_ts = add_b128(total_sum_ts, shift_bytes_right_b128<4>(total_sum_ts));
  total_sum_ts = add_b128(total_sum_ts, shift_bytes_right_b128<8>(total_sum_ts));
  total_sum_pyrts = add_b128(total_sum_pyrts, shift_bytes_right_b128<4>(total_sum_pyrts));
  total_sum_pyrts = add_b128(total_sum_pyrts, shift_bytes_right_b128<8>(total_sum_pyrts));
  total_sum_tv = add_b128(total_sum_tv, shift_bytes_right_b128<4>(total_sum_tv));
  total_sum_tv = add_b128(total_sum_tv, shift_bytes_right_b128<8>(total_sum_tv));
  total_sum_del = add_b128(total_sum_del, shift_bytes_right_b128<4>(total_sum_del));
  total_sum_del = add_b128(total_sum_del, shift_bytes_right_b128<8>(total_sum_del));
  
  TN_string_distance d= {get_int_0_b128(total_sum_del),
                         static_cast<float>(get_int_0_b128(total_sum_ts)-get_int_0_b128(total_sum_pyrts)),
                         static_cast<float>(get_int_0_b128(total_sum_pyrts)),
			 static_cast<float>(get_int_0_b128(total_sum_tv))};
  
  
  return d;
}
//---------------------------------
// Tamura Nei.
//
// For this distance function we need to count the number of
// transversions and the two types of transitions. i.e. the number of
// transitions ts1:A<->G and ts2:T<->C. Notice that in the coding of
// the nucleotides the second bit decides if it is a pyrimidine or a
// purine. While the first bit discribes what kind of pyrimidine or
// purine it is.
//
// Thus to compute the tv we just do as before:
// tv = and_b128(shift_each32_bits_right_b128(diff,ONE),LEAST_SIGNIFCANT_BIT);
//
// to compute ts we compute:
// ts = andnot_b128(tv, and_b128(diff,LEAST_SIGNIFCANT_BIT));
//
// to compute pyrimidine transitions we do:
// arepyrimidines = shift_each32_bits_right_b128(and_b128(get_b128(ptr1),get_b128(ptr2)),ONE);
// pyrts = and_b128(arepyrimidines,ts);
//
// the number of purine transitions is now given by the number of transitions - the number of pyrimidine transitions
//--------------------------------

//--------------------------------------
// LEVEL 1
//
// BLOCK SIZE 2 BITS REPRESENTING ATMOST <= 3
//
// This function reads three b128s from ptr1 and ptr2 and returns a pair
// that describes how many transitions and transversions there are in
// the three read b128s.

static void
dist_level_1(b128 &sum_ts_l1, b128 &sum_pyrts_l1, b128 &sum_tv_l1, b128 &sum_del_l1){
  b128 diff1;
  b128 diff2;
  b128 diff3;
  b128 ptrand1;
  b128 ptrand2;
  b128 ptrand3;
  b128 del1;
  b128 del2;
  b128 del3;
  b128 *_del_ptr1;
  b128 *_del_ptr2;
  b128 *_ptr1;
  b128 *_ptr2;
  
  //--------
  // LOOP 1
  // del has '11' in the blocks which should be disregarded.
  // diff = ( ~del) & (ptr1 ^ ptr2 )
  // diff has ones in the blocks that differ.
  
  _del_ptr1 = del_ptr1;
  _del_ptr2 = del_ptr2;

  del1 = or_b128(get_b128(_del_ptr1),get_b128(_del_ptr2));
  del2 = or_b128(get_b128(_del_ptr1+1),get_b128(_del_ptr2+1));
  del3 = or_b128(get_b128(_del_ptr1+2),get_b128(_del_ptr2+2));
  del_ptr1 = _del_ptr1 + 3;
  del_ptr2 = _del_ptr2 + 3;
  
  _ptr1 = ptr1;
  _ptr2 = ptr2;
  b128 p1 = get_b128(_ptr1);
  b128 p2 = get_b128(_ptr2);
  ptrand1 = shift_each32_bits_right_b128(and_b128(p1,p2),ONE);
  diff1 = andnot_b128(del1,xor_b128(p1,p2));
  p1 = get_b128(_ptr1+1);
  p2 = get_b128(_ptr2+1);
  ptrand2 = shift_each32_bits_right_b128(and_b128(p1,p2),ONE);
  diff2 = andnot_b128(del2,xor_b128(p1,p2));
  p1 = get_b128(_ptr1+2);
  p2 = get_b128(_ptr2+2);
  ptrand3 = shift_each32_bits_right_b128(and_b128(p1,p2),ONE);
  diff3 = andnot_b128(del3,xor_b128(p1,p2));


  ptr1 = _ptr1 + 3;
  ptr2 = _ptr2 + 3;
  

  // in diff for each 2 bit block
  // XY  - the bit positions in the block
  // 00  - if equal
  // 01  - if transition
  // 11 or 10 - if transversion
  //
  // Thus a transversion has occured if there is a '1' in bit X
  // and a transversion has occured if the is a '1' in bit Y but not bit X.
  // I.e.:
  // TV = X            - which is given by shifting one bit right and masking out the least significant bit.
  // TS = (~X) & Y     - Y is given by masking out the least significant bit
  //
  // the computation of pyrimidines is described above.
  //
  // The number of deletions are the number of ones given by the del mask.
  // Modernization Phase 0 (modernization_plan.md): 'register' was a
  // no-op hint ignored by all modern compilers well before this; ISO
  // C++17 removed the keyword entirely, making it a hard error. Removed
  // here as a required syntax fix for the whole project to compile
  // under C++17 - DNA_b128 itself is out of scope for style
  // modernization this phase (see modernization_plan.md).
  b128 _sum_del_l1;
  b128 _sum_tv_l1;
  b128 _sum_ts_l1;
  b128 _sum_pyrts_l1;
  b128 tmp_tv;
  b128 tmp_ts;
  _sum_tv_l1 = and_b128(shift_each32_bits_right_b128(diff1,ONE),LEAST_SIGNIFCANT_BIT);
  _sum_del_l1 = and_b128(del1,LEAST_SIGNIFCANT_BIT);
  _sum_ts_l1 = andnot_b128(_sum_tv_l1, and_b128(diff1,LEAST_SIGNIFCANT_BIT));

  _sum_pyrts_l1 = and_b128(ptrand1,_sum_ts_l1);
  
  //
  tmp_tv = and_b128(shift_each32_bits_right_b128(diff2,ONE),LEAST_SIGNIFCANT_BIT);
  _sum_del_l1 = add_b128(_sum_del_l1, and_b128(del2,LEAST_SIGNIFCANT_BIT));
  tmp_ts = andnot_b128(tmp_tv, and_b128(diff2,LEAST_SIGNIFCANT_BIT)); 
  _sum_tv_l1 = add_b128(_sum_tv_l1, tmp_tv);
  _sum_ts_l1 = add_b128(_sum_ts_l1,tmp_ts);

  _sum_pyrts_l1 = add_b128(_sum_pyrts_l1,and_b128(ptrand2,tmp_ts));

  //
  tmp_tv = and_b128(shift_each32_bits_right_b128(diff3,ONE),LEAST_SIGNIFCANT_BIT);
  _sum_del_l1 = add_b128(_sum_del_l1, and_b128(del3,LEAST_SIGNIFCANT_BIT));
  tmp_ts = andnot_b128(tmp_tv, and_b128(diff3,LEAST_SIGNIFCANT_BIT)); 
  _sum_tv_l1 = add_b128(_sum_tv_l1, tmp_tv);
  _sum_ts_l1 = add_b128(_sum_ts_l1,tmp_ts);

  _sum_pyrts_l1 = add_b128(_sum_pyrts_l1,and_b128(ptrand3,tmp_ts));


  
  //set the locals to the globals
  sum_del_l1 = _sum_del_l1;
  sum_tv_l1 = _sum_tv_l1;
  sum_ts_l1 = _sum_ts_l1;
  sum_pyrts_l1 = _sum_pyrts_l1;

}



//---------------------------------
//LEVEL 2
//
// BLOCK SIZE 4 BITS REPRESENTING ATMOST <= 15
//
void
dist_level_2(b128 &sum_ts_l2, b128 &sum_pyrts_l2, b128 &sum_tv_l2, b128 &sum_del_l2){
  b128 sum_ts_l1;
  b128 sum_pyrts_l1;
  b128 sum_tv_l1;
  b128 sum_del_l1;
  PREFETCH_DATA(6);
  //----
  // LOOP 1
  // Call level 1 and add the blocks of size 2 into
  // a block of size 4.
  dist_level_1(sum_ts_l1, sum_pyrts_l1, sum_tv_l1, sum_del_l1);

  convert_sum(sum_ts_l2,sum_ts_l1, TWO_BIT_MASK, TWO);
  convert_sum(sum_pyrts_l2,sum_pyrts_l1, TWO_BIT_MASK, TWO);
  convert_sum(sum_tv_l2,sum_tv_l1, TWO_BIT_MASK, TWO);
  convert_sum(sum_del_l2,sum_del_l1, TWO_BIT_MASK, TWO);
  
  //----
  // LOOP 2
  dist_level_1(sum_ts_l1, sum_pyrts_l1, sum_tv_l1, sum_del_l1);
  sum_with_previous_level(sum_ts_l2,sum_ts_l1, TWO_BIT_MASK, TWO);
  sum_with_previous_level(sum_pyrts_l2,sum_pyrts_l1, TWO_BIT_MASK, TWO);
  sum_with_previous_level(sum_tv_l2,sum_tv_l1, TWO_BIT_MASK, TWO);
  sum_with_previous_level(sum_del_l2,sum_del_l1, TWO_BIT_MASK, TWO);
  
}


//---------------------------------
//LEVEL 3
//
// BLOCK SIZE 8 BITS REPRESENTING ATMOST <= 255
//
// After this function has been called 3*5*17=255 b128s have been compaired
// and the number of TS and TV are in the associated variables for
// level 3.
void
dist_level_3(b128 &sum_ts_l3, b128 &sum_pyrts_l3, b128 &sum_tv_l3, b128 &sum_del_l3){
  b128 sum_ts_l2;
  b128 sum_pyrts_l2;
  b128 sum_tv_l2;
  b128 sum_del_l2;
  //----
  // LOOP 1
  // Call level 2 and add the blocks of size 4 into
  // a block of size 8.
  dist_level_2(sum_ts_l2, sum_pyrts_l2, sum_tv_l2, sum_del_l2);
  convert_sum(sum_ts_l3,sum_ts_l2, FOUR_BIT_MASK, FOUR);
  convert_sum(sum_pyrts_l3,sum_pyrts_l2, FOUR_BIT_MASK, FOUR);
  convert_sum(sum_tv_l3,sum_tv_l2, FOUR_BIT_MASK, FOUR);
  convert_sum(sum_del_l3,sum_del_l2, FOUR_BIT_MASK, FOUR);
  

  //----
  // LOOP 2
  dist_level_2(sum_ts_l2, sum_pyrts_l2, sum_tv_l2, sum_del_l2);
  sum_with_previous_level(sum_ts_l3,sum_ts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_pyrts_l3,sum_pyrts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_tv_l3,sum_tv_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_del_l3,sum_del_l2, FOUR_BIT_MASK, FOUR);
  

  
  //----
  // LOOP 3
  dist_level_2(sum_ts_l2, sum_pyrts_l2, sum_tv_l2, sum_del_l2);
  sum_with_previous_level(sum_ts_l3,sum_ts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_pyrts_l3,sum_pyrts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_tv_l3,sum_tv_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_del_l3,sum_del_l2, FOUR_BIT_MASK, FOUR);
  
  
  //----
  // LOOP 4
  dist_level_2(sum_ts_l2, sum_pyrts_l2, sum_tv_l2, sum_del_l2);
  sum_with_previous_level(sum_ts_l3,sum_ts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_pyrts_l3,sum_pyrts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_tv_l3,sum_tv_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_del_l3,sum_del_l2, FOUR_BIT_MASK, FOUR);
  

  
  //----
  // LOOP 5
  dist_level_2(sum_ts_l2, sum_pyrts_l2, sum_tv_l2, sum_del_l2);
  sum_with_previous_level(sum_ts_l3,sum_ts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_pyrts_l3,sum_pyrts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_tv_l3,sum_tv_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_del_l3,sum_del_l2, FOUR_BIT_MASK, FOUR);
  
  
  //----
  // LOOP 6
  dist_level_2(sum_ts_l2, sum_pyrts_l2, sum_tv_l2, sum_del_l2);
  sum_with_previous_level(sum_ts_l3,sum_ts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_pyrts_l3,sum_pyrts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_tv_l3,sum_tv_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_del_l3,sum_del_l2, FOUR_BIT_MASK, FOUR);
  
  
  //----
  // LOOP 7
  dist_level_2(sum_ts_l2, sum_pyrts_l2, sum_tv_l2, sum_del_l2);
  sum_with_previous_level(sum_ts_l3,sum_ts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_pyrts_l3,sum_pyrts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_tv_l3,sum_tv_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_del_l3,sum_del_l2, FOUR_BIT_MASK, FOUR);
  
  
  //----
  // LOOP 8
  dist_level_2(sum_ts_l2, sum_pyrts_l2, sum_tv_l2, sum_del_l2);
  sum_with_previous_level(sum_ts_l3,sum_ts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_pyrts_l3,sum_pyrts_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_tv_l3,sum_tv_l2, FOUR_BIT_MASK, FOUR);
  sum_with_previous_level(sum_del_l3,sum_del_l2, FOUR_BIT_MASK, FOUR);
  
  
}//END LEVEL 3





//---------------------------------
//LEVEL 3
//
// BLOCK SIZE 16 BITS REPRESENTING ATMOST <= 65,535
//
//
void
dist_level_4(b128 &sum_ts_l4, b128 &sum_pyrts_l4, b128 &sum_tv_l4, b128 &sum_del_l4){
  b128 sum_ts_l3;
  b128 sum_pyrts_l3;
  b128 sum_tv_l3;
  b128 sum_del_l3;
  //----
  // LOOP 1
  // Call level 3 and add the blocks of size 8 into
  // a block of size 16.
  dist_level_3(sum_ts_l3, sum_pyrts_l3, sum_tv_l3, sum_del_l3);
  //  convert_sum(sum_ts_l4,sum_ts_l3, EIGHT_BIT_MASK, EIGHT);
  //  convert_sum(sum_tv_l4,sum_tv_l3, EIGHT_BIT_MASK, EIGHT);
  //  convert_sum(sum_del_l4,sum_del_l3, EIGHT_BIT_MASK, EIGHT);
  CONVERT_SUM_IMMEDIATEBYTESHIFT(sum_ts_l4,sum_ts_l3, EIGHT_BIT_MASK, 1);
  CONVERT_SUM_IMMEDIATEBYTESHIFT(sum_pyrts_l4,sum_pyrts_l3, EIGHT_BIT_MASK, 1);
  CONVERT_SUM_IMMEDIATEBYTESHIFT(sum_tv_l4,sum_tv_l3, EIGHT_BIT_MASK, 1);
  CONVERT_SUM_IMMEDIATEBYTESHIFT(sum_del_l4,sum_del_l3, EIGHT_BIT_MASK, 1);
  

  //----
  // LOOP 2-128
  //PENDING can of course unroll this loop to
  for ( int i = 127 ; i != 0 ; i-- ){
    dist_level_3(sum_ts_l3, sum_pyrts_l3, sum_tv_l3, sum_del_l3);
    //    sum_with_previous_level(sum_ts_l4,sum_ts_l3, EIGHT_BIT_MASK, EIGHT);
    //    sum_with_previous_level(sum_tv_l4,sum_tv_l3, EIGHT_BIT_MASK, EIGHT);
    //    sum_with_previous_level(sum_del_l4,sum_del_l3, EIGHT_BIT_MASK, EIGHT);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(sum_ts_l4,sum_ts_l3, EIGHT_BIT_MASK, 1);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(sum_pyrts_l4,sum_pyrts_l3, EIGHT_BIT_MASK, 1);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(sum_tv_l4,sum_tv_l3, EIGHT_BIT_MASK, 1);
    SUM_WITH_PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT(sum_del_l4,sum_del_l3, EIGHT_BIT_MASK, 1);
    
  }
}










