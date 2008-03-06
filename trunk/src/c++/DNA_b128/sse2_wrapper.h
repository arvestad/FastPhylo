//--------------------------------------------------
//                                        
// File: sse2_wrapper.h
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: sse2_wrapper.h,v 1.3 2006/12/20 09:58:19 isaac Exp $                                 
//
//--------------------------------------------------


//
// This file contains wrappers to the 128 SSE2 integer registers.
//

#ifndef SSE2_WRAPPER_H
#define SSE2_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

#include <emmintrin.h>


typedef __m128i b128;

typedef union {
  int _mem[4];
  b128 _b128;
} union_b128_mem;

//-------------------------
// GET AND SET OPERATIONS
// The returned b128 is {i3,i2,i1,i0}
static __inline b128
set_ints_b128(int i3, int i2, int i1, int i0){
  return _mm_set_epi32(i3,i2,i1,i0);
}

// The return b128 is {0,0,0,i0}
static __inline b128
set_first_int_b128(int i0){
  return _mm_cvtsi32_si128(i0);
}

static __inline b128
set_zero_b128(){
  return _mm_setzero_si128();
}

//sets all bytes to be b.
static __inline b128
set_all_bytes(char b){
  return _mm_set1_epi8(b);
}
//sets all shorts to be b
static __inline b128
set_all_shorts(short b){
  return _mm_set1_epi16(b);
}
//sets all ints to be b
static __inline b128
set_all_ints(int b){
  return _mm_set1_epi32(b);
}

//---------------------------------
// IMMEDIATE ACCESS
// Note! The pos has to be an immediate, i.e., a constant integer.
#define get_immediate_int16_b128(__A, __IMM) _mm_extract_epi16(__A,__IMM)
#define set_immediate_int16_b128(__A, __B, __IMM) _mm_insert_epi16(__A,__B,__IMM)


//------------------
//IMMEDIATE GETS
static __inline int
get_int_0_b128(b128 a){
  return  ((_mm_extract_epi16(a,1)<<16) |  _mm_extract_epi16(a,0));
}

static __inline int
get_int_1_b128(b128 a){
  return  ((_mm_extract_epi16(a,3)<<16) |  _mm_extract_epi16(a,2));
}

static __inline int
get_int_2_b128(b128 a){
  return  ((_mm_extract_epi16(a,5)<<16) |  _mm_extract_epi16(a,4));
}

static __inline int
get_int_3_b128(b128 a){
  return  ((_mm_extract_epi16(a,7)<<16) |  _mm_extract_epi16(a,6));
}
//------------------
// IMMEDIATE SETS
static __inline b128
set_int_0_b128(b128 a, int b){
  return _mm_insert_epi16(_mm_insert_epi16(a, (b>>16) ,1),b,0);
}

static __inline b128
set_int_1_b128(b128 a, int b){
  return _mm_insert_epi16(_mm_insert_epi16(a, (b>>16) ,3),b,2);
}

static __inline b128
set_int_2_b128(b128 a, int b){
  return _mm_insert_epi16(_mm_insert_epi16(a, (b>>16) ,5),b,4);
}

static __inline b128
set_int_3_b128(b128 a, int b){
  return _mm_insert_epi16(_mm_insert_epi16(a, (b>>16) ,7),b,6);
}

//----------------------------------
// Returns the 32 bits at position pos in the 128 bits.
// E.g. get_int_b128(set_b128(0,0,X,0),1) returns X.
static __inline int
get_int_b128(b128 a, int intpos){
  union_b128_mem b;
  b._b128 = a;
  return b._mem[intpos];
}

static __inline b128
set_int_b128(b128 a, int b, int intpos){
  union_b128_mem c;
  c._b128 = a;
  c._mem[intpos] = b;

  return c._b128;
}


//--------------------------------
// GETTING AND SETTING OF BITS
static __inline int
get_bit_b128(b128 a, int bitpos){
  // bitpos is in the interval [0,128) and consists of 7 bits:
  //the first five bits is an index of an int
  //the last two is one of the four ints.
  return (( get_int_b128(a,bitpos>>5) >> (bitpos & 0x1F)) & 0x1);
}

static __inline b128
set_bit_b128(b128 a, int bitpos, int bitvalue){
  union_b128_mem b;
  
  b._b128 = a;
  b._mem[bitpos>>5] |= (0x1 <<  (bitpos & 0x1F) );

  return b._b128;
}




//-------------------------
// COMPAIR
// Compair if all bits are equal.
// If not equal 0 is returned if equal 0xFFffFFff is returned
static __inline int
equal_b128(b128 a, b128 b){
  b128 c = _mm_cmpeq_epi32(a,b);

  return get_int_0_b128(c) | get_int_1_b128(c) |
    get_int_2_b128(c) | get_int_3_b128(c);
}

//--------------------------
// ADDITION
static __inline b128
add_b128(b128 a, b128 b){
  return _mm_add_epi32(a,b);
}
// SUBTRACTION
static __inline b128
sub_b128(b128 a, b128 b){
  return _mm_sub_epi32(a,b);
}

// BITWISE XOR returns a ^ b
static __inline b128
xor_b128(b128 a, b128 b){
  return _mm_xor_si128(a,b);
}
// BITWISE OR returns a | b
static __inline b128
or_b128(b128 a, b128 b){
  return _mm_or_si128(a,b);
}
// BITWISE AND returns a & b
static __inline b128
and_b128(b128 a, b128 b){
  return _mm_and_si128(a,b);
}
// BITWISE AND NOT returns (~a) & b
static __inline b128
andnot_b128(b128 a, b128 b){
  return _mm_andnot_si128(a,b);
}
//BITWISE NEGATE returns (~a)
static __inline b128
negate_b128(b128 a){
  return xor_b128(_mm_set1_epi8(0xff),a);
}

//-------------------------
// SHIFTING ALL 128 BITS IMMEDIATE BYTES
// NOTEsse2_wrapper! An immediate is a constant integer (literal),
// thus immediate_bytes can not be a variable argument!
// Shift all 128 bits in a by specific number of bytes
// E.g. get_int(
//              shift_bytes_left_b128(set_b128(0,0,X,0),4)
//             ,2) returns X.
#define shift_bytes_left_b128(__A, __IMM_BYTES)  _mm_slli_si128(__A,__IMM_BYTES)

//E.g. get_int(
//             shift_bytes_right_b128(set_b128(0,0,X,0),4)
//            ,0) returns X.
#define shift_bytes_right_b128(__A, __IMM_BYTES)  _mm_srli_si128(__A,__IMM_BYTES)

//-------------------------
// LOCAL SHIFTING
// Local shifting of each 32 bits zeros are added to the end.
// E.g. get_int(
//             shift_each32_bits_right_b128(set_b128(0,0,2,0),set_first_int_b128(1))
//             ,1) returns 1
static __inline b128
shift_each32_bits_right_b128(b128 a, b128 b){
  return _mm_srl_epi32(a,b);
}

// E.g. get_int(
//             shift_each32_bits_left_b128(set_b128(0,0,2,0),set_first_int_b128(1))
//             ,1) returns 4
static __inline b128
shift_each32_bits_left_b128(b128 a, b128 b){
  return _mm_sll_epi32(a,b);
}
     


//---------------------------
// MEMORY ALLOC/FREE
// Allocates num_b128*16 bytes memory that is 16 bytes aligned.

b128 *alloc_b128(int num_b128);
b128 *calloc_b128(int num_b128);

void free_b128(b128 *a);

//-----------------
// MEMORY ACCESS
// Note! Accessing the memory as ptr[pos] works but
// I am not sure wether the memory is assumed to be aligned.
// Therefore maybe *_loadu_* is used instead, which of course
// is slower than *_load_*.
static __inline b128
get_from_mem_b128(const b128 *ptr,int pos){
  return _mm_load_si128(ptr+pos);
}
static __inline b128
get_b128(const b128 *ptr){
  return _mm_load_si128(ptr);
}

static __inline void
set_b128(b128 *ptr, b128 val){
  _mm_store_si128(ptr, val);
}

static __inline void
set_b128_NOCACHEPOLLUTE(b128 *ptr, b128 val){
  _mm_stream_si128(ptr, val);
}


//-----------------------------------------------------
// getting and setting of bits in aligned memory
// The bitpos can be any in the interval [0,128*num_b128)
static __inline int
get_bit_in_mem_b128(const b128 *ptr, int bitpos){
  //take the least significant seven bits and send them to the get_bit function
  return get_bit_b128(get_from_mem_b128(ptr,(bitpos>>7)),
                      (bitpos & 0x7F) );
}

static __inline void
set_bit_in_mem_b128(b128 *ptr, int bitpos, int bitvalue){
  //take the least significant seven bits and send them to the get_bit function
  b128 *b128pos = ptr+(bitpos>>7);
  _mm_store_si128(b128pos,
                  set_bit_b128(_mm_load_si128(b128pos),
                               (bitpos & 0x7F),
                               bitvalue ));
}


//----------------------------------
//LOGICAL OPERATIONS ON ARRAYS
//
// dst <- dst OR src
static __inline void
mem_or_b128(b128 *dst, const b128 *src, const int num_b128s){

    b128 t1, t2;
    //walk through four 128s in each row
    const b128 *src_end = src+(num_b128s-3);

    while ( src < src_end ){
      _mm_prefetch((const char*)(src)+512,  _MM_HINT_NTA);
      _mm_prefetch((const char*)(dst)+512,  _MM_HINT_T0);
    
        t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = or_b128(t1, t2);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = or_b128(t1, t2);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = or_b128(t1, t2);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = or_b128(t1, t2);
        set_b128(dst++, t1);    
    } 
    
    //compute or of remaining b128s,
    while ( src < src_end+3){
        t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = or_b128(t1, t2);
        set_b128(dst++, t1);
    }
}

static __inline void
mem_and_b128(b128 *dst, const b128 *src, const int num_b128s){
    b128 t1, t2;
    //walk through four 128s in each row
    const b128 *src_end = src+(num_b128s-3);

    while ( src < src_end ){
        _mm_prefetch((const char*)(src)+512,  _MM_HINT_NTA);
	_mm_prefetch((const char*)(dst)+512,  _MM_HINT_T0);
    
        t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = and_b128(t1, t2);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = and_b128(t1, t2);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = and_b128(t1, t2);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = and_b128(t1, t2);
        set_b128(dst++, t1);    
    } 
    
    //compute xor of remaining b128s,
    while ( src < src_end+3){
        t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = and_b128(t1, t2);
        set_b128(dst++, t1);
    }
}

static __inline void
mem_xor_b128(b128 *dst, const b128 *src, const int num_b128s){
    b128 t1, t2;
    //walk through four 128s in each row
    const b128 *src_end = src+(num_b128s-3);

    while ( src < src_end ){
        _mm_prefetch((const char*)(src)+512,  _MM_HINT_NTA);
	_mm_prefetch((const char*)(dst)+512,  _MM_HINT_T0);
    
        t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = xor_b128(t1, t2);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = xor_b128(t1, t2);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = xor_b128(t1, t2);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = xor_b128(t1, t2);
        set_b128(dst++, t1);    
    } 
    
    //compute xor of remaining b128s,
    while ( src < src_end+3){
        t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = xor_b128(t1, t2);
        set_b128(dst++, t1);
    }
}

// dst <- (~dst) & src
static __inline void
mem_andnot_b128(b128 *dst, const b128 *src, const int num_b128s){
    b128 t1, t2;
    //walk through four 128s in each row
    const b128 *src_end = src+(num_b128s-3);

    while ( src < src_end ){
        _mm_prefetch((const char*)(src)+512,  _MM_HINT_NTA);
	_mm_prefetch((const char*)(dst)+512,  _MM_HINT_T0);
    
        t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = andnot_b128(t1, t2);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = andnot_b128(t1, t2);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = andnot_b128(t1, t2);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = andnot_b128(t1, t2);
        set_b128(dst++, t1);    
    } 
    
    //compute xor of remaining b128s,
    while ( src < src_end+3){
        t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = andnot_b128(t1, t2);
        set_b128(dst++, t1);
    }
}

//dst <- dst & (~src)
static __inline void
mem_reversed_andnot_b128(b128 *dst, const b128 *src, const int num_b128s){

    b128 t1, t2;
    //walk through four 128s in each row
    const b128 *src_end = src+(num_b128s-3);

    while ( src < src_end ){
        _mm_prefetch((const char*)(src)+512,  _MM_HINT_NTA);
	_mm_prefetch((const char*)(dst)+512,  _MM_HINT_T0);    
        t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = andnot_b128(t2, t1);//reversed
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = andnot_b128(t2, t1);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = andnot_b128(t2, t1);
        set_b128(dst++, t1);
        
	t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = andnot_b128(t2, t1);
        set_b128(dst++, t1);    
    } 
    
    //compute xor of remaining b128s,
    while ( src < src_end+3){
        t1 = get_b128(src++);
        t2 = get_b128(dst);
        t1 = andnot_b128(t2, t1);
        set_b128(dst++, t1);
    }
}
//---------------------------
// PRINTING
void print_ints_b128(b128 a);
void print_bits_32(int a);
void print_bits_b128(b128 a);

//Divides a into blocks of size block_size and prints the integer in
//each such block. block_size has to be a power of 2.
void print_blocks_b128(b128 a, int block_size);


  // PENTIUM SPECIFIC TICKS COUNTER
  typedef unsigned long long ticks;

  static __inline__ ticks getticks(void)
  {
    unsigned a, d; 
    asm volatile("rdtsc" : "=a" (a), "=d" (d)); 
    return ((ticks)a) | (((ticks)d) << 32); 
  }


  
#ifdef __cplusplus
}
#endif  

#endif // SSE2_WRAPPER_H
