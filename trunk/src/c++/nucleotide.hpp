//--------------------------------------------------
//                                        
// File: nucleotide.hpp
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: nucleotide.hpp,v 1.6 2006/12/27 14:20:49 isaac Exp $                                 
//
//--------------------------------------------------


#ifndef NUCLEOTIDE_HPP
#define NUCLEOTIDE_HPP


#include <assert.h>

//-------------------------------------------------------------
// TERMS  GLOSSARY
//
// DNA: Adenine, Cytosine, Guanine, Thymidine
//  A pairs with T
//  C pairs with G
//
// Purine: An Adenine (A) or a Guanine (G) nucleotide.
//
// Pyrimidine: A Cytosine (C) or a Thymidine (T) nucleotide. 
//
// RNA: A, C, G, Uracil
//  DNA side   RNA side
//  A            U
//  T            A
//  C            G
//  G            C
//
// Transition: (Ti) This is a substitution of a purine for a purine or
// a pyrimidine for a pyrimidine, i.e. A for G and C for T.
//
// Transversion: (Tv) This is a substitution of a purine for a
// pyrimidine or vice versa.
//
// Transition-transversion ratio (Ti/Tv): Is often 2, i.e., transitions occur
// more frequently than transversions.
//
//-------------------------------------------------------------------


enum nucleotide {
  //IUPAC Codes and phylip codes
  DNA_A_ = 0, //'A', 'a'
  DNA_C_ = 3, //'C', 'c'
  DNA_G_ = 1, //'G', 'g'
  DNA_T_ = 2, //'T', 't', 'U', 'u'  
  //AMBIGUITY CODES           COMPLEMENT CODE
  DNA_UNKNOWN_ = 4, //'.', or '-' (these occur in alignments)
  DNA_NOT_ALLOWED = 5,
  // The 6th bit is a flag signifying that it is an ambiguity
  // The four bits in positions 0,1,2,3 define which the possible symbols are (see ambiguity flags below).
  // 1 in pos 0 means that A is possible
  // 1 in pos 1 means that C is possible
  // 1 in pos 2 means that G is possible
  // 1 in pos 3 means that T is possible
  DNA_M_ = 6,// 	A or C 	        K
  DNA_R_ = 7,// 	A or G 	        Y
  DNA_W_ = 8,// 	A or T 	        W
  DNA_S_ = 9,// 	C or G 	        S
  DNA_Y_ = 10,// 	C or T 	        R
  DNA_K_ = 11,// 	G or T 	        M
  DNA_V_ = 12,// 	A or C or G 	B
  DNA_H_ = 13,// 	A or C or T 	D
  DNA_D_ = 14,// 	A or G or T 	H
  DNA_B_ = 15,// 	C or G or T 	V
  DNA_N_ = 16,//  A or C or G or T  // 'N', 'X', 'x', or '?',

//   // Ambiguity flags (use bitwise or of these to create ambiguities)
//   // These are not real nucleotides! ONLY BITWISE FLAGS! 
//   AMBIG_A_FLAG = 0x21,
//   AMBIG_C_FLAG = 0x22,
//   AMBIG_G_FLAG = 0x24,
//   AMBIG_T_FLAG = 0x28,
};


//Note that xor for a transition is 0x01, i.e.,
// A xor G = 0x01
// C xor T = 0x01
//and that xor for a transversion are either 0x02 or 0x03.


#define NUCLEOTIDE_INT_MASK 0x3

//returns !=0 if the nucleotide is A,C,G, or T
//and 0 if the nucletide is something else
static __inline int
is_regular_nucleotide(const nucleotide &n){
  return (n < 4);
}

static __inline int
is_ambiguity_nucleotide(const nucleotide &n){
  return (n > 5);
}



//-----------------------------
// NUCLEOTIDE DIFF
// Notice that in the coding of the nucleotides the second bit decides
// if it is a pyrimidine or a purine. While the first bit discribes
// what kind of pyrimidine or purine it is.

enum nucleotide_difference{
  EQUAL, //if a,g,c,t equal
  PURINE_TRANSITION,
  PYRIMIDINE_TRANSITION,
  TRANSVERSION,
  DELETION,//if one is an unknown
  AMBIGIOUS//if one is ambigious
};

static __inline nucleotide_difference
get_nucleotide_difference(nucleotide n1, nucleotide n2){
  if ( is_regular_nucleotide(n1) && is_regular_nucleotide(n2) ){//if both regular
    int tmp = n1 ^ n2;
    if ( tmp ==  0 )
      return EQUAL;
    if ( tmp ==  1 ) {// transition
      if( (n1 & 0x2) == 0 )
        return PURINE_TRANSITION;
      else
        return PYRIMIDINE_TRANSITION;
    }
    else
      return TRANSVERSION;
  }
  
  if ( n1==DNA_UNKNOWN_ || n2==DNA_UNKNOWN_)  
    return DELETION;

  assert( is_ambiguity_nucleotide(n1) || is_ambiguity_nucleotide(n2) );

  return AMBIGIOUS;
}
//---------------

//-------------------------------------------------------
// LOOK UP TABLE

const static nucleotide CHAR2NUCLEOTIDE_MAP[256] = {
  //---------------------------
  DNA_NOT_ALLOWED /* 0x0 */ ,DNA_NOT_ALLOWED /* 0x1 */ ,DNA_NOT_ALLOWED /* 0x2 */ ,DNA_NOT_ALLOWED /* 0x3 */ ,DNA_NOT_ALLOWED /* 0x4 */ ,DNA_NOT_ALLOWED /* 0x5 */ ,DNA_NOT_ALLOWED /* 0x6 */ ,DNA_NOT_ALLOWED /* 0x7 */ ,DNA_NOT_ALLOWED /* 0x8 */ ,DNA_NOT_ALLOWED /* 0x9 */ ,DNA_NOT_ALLOWED /* 0xa */ ,DNA_NOT_ALLOWED /* 0xb */ ,DNA_NOT_ALLOWED /* 0xc */ ,DNA_NOT_ALLOWED /* 0xd */ ,DNA_NOT_ALLOWED /* 0xe */ ,DNA_NOT_ALLOWED /* 0xf */ ,DNA_NOT_ALLOWED /* 0x10 */ ,DNA_NOT_ALLOWED /* 0x11 */ ,DNA_NOT_ALLOWED /* 0x12 */ ,DNA_NOT_ALLOWED /* 0x13 */ ,DNA_NOT_ALLOWED /* 0x14 */ ,DNA_NOT_ALLOWED /* 0x15 */ ,DNA_NOT_ALLOWED /* 0x16 */ ,DNA_NOT_ALLOWED /* 0x17 */ ,DNA_NOT_ALLOWED /* 0x18 */ ,DNA_NOT_ALLOWED /* 0x19 */ ,DNA_NOT_ALLOWED /* 0x1a */ ,DNA_NOT_ALLOWED /* 0x1b */ ,DNA_NOT_ALLOWED /* 0x1c */ ,DNA_NOT_ALLOWED /* 0x1d */ ,DNA_NOT_ALLOWED /* 0x1e */ ,DNA_NOT_ALLOWED /* 0x1f */ ,DNA_NOT_ALLOWED /* 0x20 */ ,DNA_NOT_ALLOWED /* 0x21 */ ,DNA_NOT_ALLOWED /* 0x22 */ ,DNA_NOT_ALLOWED /* 0x23 */ ,DNA_NOT_ALLOWED /* 0x24 */ ,DNA_NOT_ALLOWED /* 0x25 */ ,DNA_NOT_ALLOWED /* 0x26 */ ,DNA_NOT_ALLOWED /* 0x27 */ ,DNA_NOT_ALLOWED /* 0x28 */ ,DNA_NOT_ALLOWED /* 0x29 */ ,DNA_NOT_ALLOWED /* 0x2a */ ,DNA_NOT_ALLOWED /* 0x2b */ ,DNA_NOT_ALLOWED /* 0x2c */ ,
  //-------------------------------
  DNA_UNKNOWN_ /* 0x2d  '-' */ ,
  DNA_UNKNOWN_ /* 0x2e  '.' */ ,
  //-------------------------------
  DNA_NOT_ALLOWED /* 0x2f */ ,DNA_NOT_ALLOWED /* 0x30  '0' */ ,DNA_NOT_ALLOWED /* 0x31  '1' */ ,DNA_NOT_ALLOWED /* 0x32  '2' */ ,DNA_NOT_ALLOWED /* 0x33  '3' */ ,DNA_NOT_ALLOWED /* 0x34  '4' */ ,DNA_NOT_ALLOWED /* 0x35  '5' */ ,DNA_NOT_ALLOWED /* 0x36  '6' */ ,DNA_NOT_ALLOWED /* 0x37  '7' */ ,DNA_NOT_ALLOWED /* 0x38  '8' */ ,DNA_NOT_ALLOWED /* 0x39  '9' */ ,DNA_NOT_ALLOWED /* 0x3a */ ,DNA_NOT_ALLOWED /* 0x3b */ ,DNA_NOT_ALLOWED /* 0x3c */ ,DNA_NOT_ALLOWED /* 0x3d */ ,DNA_NOT_ALLOWED /* 0x3e */ ,
  //-----------------------------
  DNA_N_ /* 0x3f  '?' */ ,
  //--------------------------
  DNA_NOT_ALLOWED /* 0x40 */ ,
  //--------------------------
  DNA_A_ /* 0x41  'A' */ ,
  DNA_B_ /* 0x42  'B' */ ,
  DNA_C_ /* 0x43  'C' */ ,
  DNA_D_ /* 0x44  'D' */ ,
  DNA_NOT_ALLOWED /* 0x45  'E' */ ,
  DNA_NOT_ALLOWED /* 0x46  'F' */ ,
  DNA_G_ /* 0x47  'G' */ ,
  DNA_H_ /* 0x48  'H' */ ,
  DNA_NOT_ALLOWED /* 0x49  'I' */ ,
  DNA_NOT_ALLOWED /* 0x4a  'J' */ ,
  DNA_K_ /* 0x4b  'K' */ ,
  DNA_NOT_ALLOWED /* 0x4c  'L' */ ,
  DNA_M_ /* 0x4d  'M' */ ,
  DNA_N_ /* 0x4e  'N' */ ,
  DNA_NOT_ALLOWED /* 0x4f  'O' */ ,
  DNA_NOT_ALLOWED /* 0x50  'P' */ ,
  DNA_NOT_ALLOWED /* 0x51  'Q' */ ,
  DNA_R_ /* 0x52  'R' */ ,
  DNA_S_ /* 0x53  'S' */ ,
  DNA_T_ /* 0x54  'T' */ ,
  DNA_T_ /* 0x55  'U' */ ,
  DNA_V_ /* 0x56  'V' */ ,
  DNA_W_ /* 0x57  'W' */ ,
  DNA_N_ /* 0x58  'X' */ ,
  DNA_Y_ /* 0x59  'Y' */ ,
  DNA_NOT_ALLOWED /* 0x5a  'Z' */ ,
  //--------------------------
  DNA_NOT_ALLOWED /* 0x5b */ ,DNA_NOT_ALLOWED /* 0x5c */ ,DNA_NOT_ALLOWED /* 0x5d */ ,DNA_NOT_ALLOWED /* 0x5e */ ,DNA_NOT_ALLOWED /* 0x5f */ ,DNA_NOT_ALLOWED /* 0x60 */ ,
  //--------------------------
  DNA_A_ /* 0x61  'a' */ ,
  DNA_B_ /* 0x62  'b' */ ,
  DNA_C_ /* 0x63  'c' */ ,
  DNA_D_ /* 0x64  'd' */ ,
  DNA_NOT_ALLOWED /* 0x65  'e' */ ,
  DNA_NOT_ALLOWED /* 0x66  'f' */ ,
  DNA_G_ /* 0x67  'g' */ ,
  DNA_H_ /* 0x68  'h' */ ,
  DNA_NOT_ALLOWED /* 0x69  'i' */ ,
  DNA_NOT_ALLOWED /* 0x6a  'j' */ ,
  DNA_K_ /* 0x6b  'k' */ ,
  DNA_NOT_ALLOWED /* 0x6c  'l' */ ,
  DNA_M_ /* 0x6d  'm' */ ,
  DNA_N_ /* 0x6e  'n' */ ,
  DNA_NOT_ALLOWED /* 0x6f  'o' */ ,
  DNA_NOT_ALLOWED /* 0x70  'p' */ ,
  DNA_NOT_ALLOWED /* 0x71  'q' */ ,
  DNA_R_ /* 0x72  'r' */ ,
  DNA_S_ /* 0x73  's' */ ,
  DNA_T_ /* 0x74  't' */ ,
  DNA_T_ /* 0x75  'u' */ ,
  DNA_V_ /* 0x76  'v' */ ,
  DNA_W_ /* 0x77  'w' */ ,
  DNA_N_ /* 0x78  'x' */ ,
  DNA_Y_ /* 0x79  'y' */ ,
  DNA_NOT_ALLOWED /* 0x7a  'z' */ ,
  //--------------------------
  DNA_NOT_ALLOWED /* 0x7b */ ,DNA_NOT_ALLOWED /* 0x7c */ ,DNA_NOT_ALLOWED /* 0x7d */ ,DNA_NOT_ALLOWED /* 0x7e */ ,DNA_NOT_ALLOWED /* 0x7f */ ,DNA_NOT_ALLOWED /* 0x80 */ ,DNA_NOT_ALLOWED /* 0x81 */ ,DNA_NOT_ALLOWED /* 0x82 */ ,DNA_NOT_ALLOWED /* 0x83 */ ,DNA_NOT_ALLOWED /* 0x84 */ ,DNA_NOT_ALLOWED /* 0x85 */ ,DNA_NOT_ALLOWED /* 0x86 */ ,DNA_NOT_ALLOWED /* 0x87 */ ,DNA_NOT_ALLOWED /* 0x88 */ ,DNA_NOT_ALLOWED /* 0x89 */ ,DNA_NOT_ALLOWED /* 0x8a */ ,DNA_NOT_ALLOWED /* 0x8b */ ,DNA_NOT_ALLOWED /* 0x8c */ ,DNA_NOT_ALLOWED /* 0x8d */ ,DNA_NOT_ALLOWED /* 0x8e */ ,DNA_NOT_ALLOWED /* 0x8f */ ,DNA_NOT_ALLOWED /* 0x90 */ ,DNA_NOT_ALLOWED /* 0x91 */ ,DNA_NOT_ALLOWED /* 0x92 */ ,DNA_NOT_ALLOWED /* 0x93 */ ,DNA_NOT_ALLOWED /* 0x94 */ ,DNA_NOT_ALLOWED /* 0x95 */ ,DNA_NOT_ALLOWED /* 0x96 */ ,DNA_NOT_ALLOWED /* 0x97 */ ,DNA_NOT_ALLOWED /* 0x98 */ ,DNA_NOT_ALLOWED /* 0x99 */ ,DNA_NOT_ALLOWED /* 0x9a */ ,DNA_NOT_ALLOWED /* 0x9b */ ,DNA_NOT_ALLOWED /* 0x9c */ ,DNA_NOT_ALLOWED /* 0x9d */ ,DNA_NOT_ALLOWED /* 0x9e */ ,DNA_NOT_ALLOWED /* 0x9f */ ,DNA_NOT_ALLOWED /* 0xa0 */ ,DNA_NOT_ALLOWED /* 0xa1 */ ,DNA_NOT_ALLOWED /* 0xa2 */ ,DNA_NOT_ALLOWED /* 0xa3 */ ,DNA_NOT_ALLOWED /* 0xa4 */ ,DNA_NOT_ALLOWED /* 0xa5 */ ,DNA_NOT_ALLOWED /* 0xa6 */ ,DNA_NOT_ALLOWED /* 0xa7 */ ,DNA_NOT_ALLOWED /* 0xa8 */ ,DNA_NOT_ALLOWED /* 0xa9 */ ,DNA_NOT_ALLOWED /* 0xaa */ ,DNA_NOT_ALLOWED /* 0xab */ ,DNA_NOT_ALLOWED /* 0xac */ ,DNA_NOT_ALLOWED /* 0xad */ ,DNA_NOT_ALLOWED /* 0xae */ ,DNA_NOT_ALLOWED /* 0xaf */ ,DNA_NOT_ALLOWED /* 0xb0 */ ,DNA_NOT_ALLOWED /* 0xb1 */ ,DNA_NOT_ALLOWED /* 0xb2 */ ,DNA_NOT_ALLOWED /* 0xb3 */ ,DNA_NOT_ALLOWED /* 0xb4 */ ,DNA_NOT_ALLOWED /* 0xb5 */ ,DNA_NOT_ALLOWED /* 0xb6 */ ,DNA_NOT_ALLOWED /* 0xb7 */ ,DNA_NOT_ALLOWED /* 0xb8 */ ,DNA_NOT_ALLOWED /* 0xb9 */ ,DNA_NOT_ALLOWED /* 0xba */ ,DNA_NOT_ALLOWED /* 0xbb */ ,DNA_NOT_ALLOWED /* 0xbc */ ,DNA_NOT_ALLOWED /* 0xbd */ ,DNA_NOT_ALLOWED /* 0xbe */ ,DNA_NOT_ALLOWED /* 0xbf */ ,DNA_NOT_ALLOWED /* 0xc0 */ ,DNA_NOT_ALLOWED /* 0xc1 */ ,DNA_NOT_ALLOWED /* 0xc2 */ ,DNA_NOT_ALLOWED /* 0xc3 */ ,DNA_NOT_ALLOWED /* 0xc4 */ ,DNA_NOT_ALLOWED /* 0xc5 */ ,DNA_NOT_ALLOWED /* 0xc6 */ ,DNA_NOT_ALLOWED /* 0xc7 */ ,DNA_NOT_ALLOWED /* 0xc8 */ ,DNA_NOT_ALLOWED /* 0xc9 */ ,DNA_NOT_ALLOWED /* 0xca */ ,DNA_NOT_ALLOWED /* 0xcb */ ,DNA_NOT_ALLOWED /* 0xcc */ ,DNA_NOT_ALLOWED /* 0xcd */ ,DNA_NOT_ALLOWED /* 0xce */ ,DNA_NOT_ALLOWED /* 0xcf */ ,DNA_NOT_ALLOWED /* 0xd0 */ ,DNA_NOT_ALLOWED /* 0xd1 */ ,DNA_NOT_ALLOWED /* 0xd2 */ ,DNA_NOT_ALLOWED /* 0xd3 */ ,DNA_NOT_ALLOWED /* 0xd4 */ ,DNA_NOT_ALLOWED /* 0xd5 */ ,DNA_NOT_ALLOWED /* 0xd6 */ ,DNA_NOT_ALLOWED /* 0xd7 */ ,DNA_NOT_ALLOWED /* 0xd8 */ ,DNA_NOT_ALLOWED /* 0xd9 */ ,DNA_NOT_ALLOWED /* 0xda */ ,DNA_NOT_ALLOWED /* 0xdb */ ,DNA_NOT_ALLOWED /* 0xdc */ ,DNA_NOT_ALLOWED /* 0xdd */ ,DNA_NOT_ALLOWED /* 0xde */ ,DNA_NOT_ALLOWED /* 0xdf */ ,DNA_NOT_ALLOWED /* 0xe0 */ ,DNA_NOT_ALLOWED /* 0xe1 */ ,DNA_NOT_ALLOWED /* 0xe2 */ ,DNA_NOT_ALLOWED /* 0xe3 */ ,DNA_NOT_ALLOWED /* 0xe4 */ ,DNA_NOT_ALLOWED /* 0xe5 */ ,DNA_NOT_ALLOWED /* 0xe6 */ ,DNA_NOT_ALLOWED /* 0xe7 */ ,DNA_NOT_ALLOWED /* 0xe8 */ ,DNA_NOT_ALLOWED /* 0xe9 */ ,DNA_NOT_ALLOWED /* 0xea */ ,DNA_NOT_ALLOWED /* 0xeb */ ,DNA_NOT_ALLOWED /* 0xec */ ,DNA_NOT_ALLOWED /* 0xed */ ,DNA_NOT_ALLOWED /* 0xee */ ,DNA_NOT_ALLOWED /* 0xef */ ,DNA_NOT_ALLOWED /* 0xf0 */ ,DNA_NOT_ALLOWED /* 0xf1 */ ,DNA_NOT_ALLOWED /* 0xf2 */ ,DNA_NOT_ALLOWED /* 0xf3 */ ,DNA_NOT_ALLOWED /* 0xf4 */ ,DNA_NOT_ALLOWED /* 0xf5 */ ,DNA_NOT_ALLOWED /* 0xf6 */ ,DNA_NOT_ALLOWED /* 0xf7 */ ,DNA_NOT_ALLOWED /* 0xf8 */ ,DNA_NOT_ALLOWED /* 0xf9 */ ,DNA_NOT_ALLOWED /* 0xfa */ ,DNA_NOT_ALLOWED /* 0xfb */ ,DNA_NOT_ALLOWED /* 0xfc */ ,DNA_NOT_ALLOWED /* 0xfd */ ,DNA_NOT_ALLOWED /* 0xfe */ ,DNA_NOT_ALLOWED /* 0xff */
};


//--------------------------------------------------
// CONVERTING 
// Uracil and Thymidine are both translated into DNA_T_
static __inline nucleotide
char2nucleotide(const char realC){
  return CHAR2NUCLEOTIDE_MAP[((int)realC)&0xff];
  
}

static __inline char
nucleotide2char(nucleotide n){
  switch (n){
  case DNA_A_: return 'a';
  case DNA_T_: return 't';
  case DNA_G_: return 'g';
  case DNA_C_: return 'c';
  case DNA_UNKNOWN_: return '-';
  case DNA_M_: return 'm';
  case DNA_R_: return 'r';
  case DNA_W_: return 'w';
  case DNA_S_: return 's';
  case DNA_Y_: return 'y';
  case DNA_K_: return 'k';
  case DNA_V_: return 'v';
  case DNA_H_: return 'h';
  case DNA_D_: return 'd';
  case DNA_B_: return 'b';
  case DNA_N_: return 'n';
  case DNA_NOT_ALLOWED: return 'Z';
  default: return 'Z';
  }
  assert( "shouldn't come here\n" == 0 );
  return 'Z';
}


#endif // NUCLEOTIDE_H
