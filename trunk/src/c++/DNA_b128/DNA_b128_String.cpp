//--------------------------------------------------
//                                        
// File: DNA_b128_String.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: DNA_b128_String.cpp,v 1.39 2006/12/28 13:17:18 isaac Exp $                                 
//
//--------------------------------------------------

#include "DNA_b128_String.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include "log_utils.hpp"
#include <algorithm>

using namespace std;

//-----------------------------------------
//  MEMORY ALLOCATION
void
DNA_b128_String::_init_mem(int capacity){

  numDatas = (128+capacity)/64;//the last allocated b128 is not used
  data = alloc_b128(numDatas);
  MEM_CHECK(data);
  unknownData = calloc_b128(numDatas);
  MEM_CHECK(unknownData);
  ambiguities.clear();

  numChars = 0;
  num_As_ = 0;
  num_Cs_ = 0;
  num_Gs_ = 0;
  num_Ts_ = 0;
  num_unknowns_ = 0;
}

void DNA_b128_String::reInitiate(size_t capacity){
  if( getTotalCapacity() < capacity ){
    _free_mem();
    _init_mem(capacity);
  }
  else{
    memset(unknownData, 0, sizeof(b128)*numDatas);//clear unknowns
    ambiguities.clear();
    numChars = 0;
    num_As_ = 0;
    num_Cs_ = 0;
    num_Gs_ = 0;
    num_Ts_ = 0;
    num_unknowns_ = 0;
  }
}

void
DNA_b128_String::_free_mem(){
  if ( data != NULL ){
    free_b128(data);
    data = NULL;
    free_b128(unknownData);
    unknownData = NULL;
  }
}

//----------------- CONSTRUCTORS -------------------------------------

DNA_b128_String::DNA_b128_String(){
  _init_mem(0);
}
DNA_b128_String::DNA_b128_String(int capacity){
  _init_mem(capacity);
}

DNA_b128_String::DNA_b128_String(const DNA_b128_String &str){
  data = NULL;
  unknownData = NULL;
  numDatas = 0;
  operator=(str);
}

void DNA_b128_String::operator=(const DNA_b128_String &str){

  if(numDatas != str.numDatas){
    _free_mem();
  
    _init_mem(str.getTotalCapacity()-1);
  }
  memcpy(data,str.data,sizeof(b128)*numDatas);
  memcpy(unknownData,str.unknownData,sizeof(b128)*numDatas);

  numChars = str.numChars;
  num_As_ = str.num_As_;
  num_Cs_ = str.num_Cs_;
  num_Gs_ = str.num_Gs_;
  num_Ts_ = str.num_Ts_;
  num_unknowns_ = str.num_unknowns_;
  ambiguities = str.ambiguities;
}

DNA_b128_String::DNA_b128_String(int capacity, const std::string &str){
  _init_mem(capacity);
  append(str);
}
DNA_b128_String::DNA_b128_String(int capacity, const char *c_str){
  _init_mem(capacity);
  append(c_str);
}


// DESTRUCTION
DNA_b128_String::~DNA_b128_String(){
  _free_mem();
}



//-----------------------------------------
// APPEND
int
DNA_b128_String::append(const char *c_str){
  size_t _num_unknowns_ = num_unknowns_;
  size_t _num_As_ = num_As_;
  size_t _num_Cs_ = num_Cs_;
  size_t _num_Gs_ = num_Gs_;
  size_t _num_Ts_ = num_Ts_;
      
  //-------------------  
  //INITIALIZE READING
  //the position of the char we are about to read
  int c_pos = getNumChars();
  //the b128 that we are working on
  b128 *c_b128_ptr = data + D_I(c_pos);
  //_mm_prefetch((char*)c_b128_ptr,_MM_HINT_T0); 
  // _mm_prefetch(c_str,_MM_HINT_NTA);  
  
  b128 c_b128;
  //if this b128 hasn't been used yet null it
  if ( first_char_in_dataslot(c_pos) )
    c_b128 = set_zero_b128();
  else
    c_b128 = get_b128(c_b128_ptr);
  //the int of the b128 we are working on
  int c_int = get_int_b128(c_b128, INT_D_I(c_pos));
  //shift character position c_pos in c_int to the least significant
  //position
  c_int = (c_int >> BIT_INT_I(c_pos) );

  //jump positions, i.e. charpositions that change b128 or int
  int next_b128_jump_pos = first_charpos_in_next_b128(c_pos);
  int next_int_jump_pos = first_charpos_in_next_int(c_pos);
  

  // read until an unallowed char is found
  for ( ; true ; c_str++ ){
    //nucleotide n = char2nucleotide(*c_str);
    //_mm_prefetch(c_str+32,_MM_HINT_NTA);  
    
    //UPDATE FREQUENCES
    switch ( *c_str ){
    case 'a': case 'A': _num_As_++; c_int = c_int | DNA_A_;break;
    case 'c': case 'C': _num_Cs_++; c_int = c_int | DNA_C_;break;
    case 'g': case 'G': _num_Gs_++; c_int = c_int | DNA_G_;break;
    case 't': case 'T': case 'u': case 'U':_num_Ts_++; c_int = c_int | DNA_T_;break;
    case ' ': continue;//skip white space
    case '-': case '.':
      _fastSetNucleotideUnknown(D_I(c_pos), INT_D_I(c_pos), BIT_INT_I(c_pos));
      //_fastSetNucleotideUnknown(c_pos);
      _num_unknowns_++;
      //n = DNA_A_;//if the position is unkown we set an adenine in its place
      c_int = c_int | DNA_A_;
      break;
    case 'm': case 'M': {
      _fastSetNucleotideUnknown(c_pos);
      ambiguities.push_back((ambiguity_nucleotide_at_position)
      {{DNA_M_, 0.5, 0.5,   0,   0},c_pos});
      c_int = c_int | DNA_A_; }
      break;// 	A or C          K
    case 'r': case 'R': {
      _fastSetNucleotideUnknown(c_pos);
      ambiguities.push_back((ambiguity_nucleotide_at_position)
      {{DNA_R_, 0.5, 0, 0.5,   0},c_pos});
      c_int = c_int | DNA_A_; }
      break;// 	A or G 	        Y
    case 'w': case 'W': {
      _fastSetNucleotideUnknown(c_pos);
      ambiguities.push_back((ambiguity_nucleotide_at_position)
      {{DNA_W_, 0.5, 0, 0, 0.5},c_pos});
      c_int = c_int | DNA_A_;  }
      break;// 	A or T 	        W
    case 's': case 'S': {
      _fastSetNucleotideUnknown(c_pos);
      ambiguities.push_back((ambiguity_nucleotide_at_position)
      {{DNA_S_, 0, 0.5, 0.5, 0},c_pos});
      c_int = c_int | DNA_A_;  }
      break;// 	C or G 	        S
    case 'y': case 'Y': {
      _fastSetNucleotideUnknown(c_pos);
      ambiguities.push_back((ambiguity_nucleotide_at_position)
      {{DNA_Y_, 0, 0.5, 0, 0.5},c_pos});
      c_int = c_int | DNA_A_; }
      break;// 	C or T 	        R
    case 'k': case 'K': {
      _fastSetNucleotideUnknown(c_pos);
      ambiguities.push_back((ambiguity_nucleotide_at_position)
      {{DNA_K_, 0, 0, 0.5, 0.5},c_pos});
      c_int = c_int | DNA_A_;  }
      break;// 	G or T 	        M
    case 'v': case 'V': {
      _fastSetNucleotideUnknown(c_pos);
      ambiguities.push_back((ambiguity_nucleotide_at_position)
      {{DNA_V_, 0.333333, 0.333333, 0.333333, 0},c_pos});
      c_int = c_int | DNA_A_;  }
      break;// 	A or C or G 	B
    case 'h': case 'H':  {
      _fastSetNucleotideUnknown(c_pos);
      ambiguities.push_back((ambiguity_nucleotide_at_position)
			    {{DNA_H_, 0.333333, 0.333333, 0, 0.333333},c_pos});
      c_int = c_int | DNA_A_;  }
      break;// 	A or C or T 	D
    case 'd': case 'D': {
      _fastSetNucleotideUnknown(c_pos);
      ambiguities.push_back((ambiguity_nucleotide_at_position)
      {{DNA_D_, 0.333333, 0, 0.333333, 0.333333},c_pos});
      c_int = c_int | DNA_A_; }
      break;// 	A or G or T 	H
    case 'b': case 'B': {
      _fastSetNucleotideUnknown(c_pos);
      ambiguities.push_back((ambiguity_nucleotide_at_position)
      {{DNA_B_, 0, 0.333333, 0.333333, 0.333333},c_pos});
      c_int = c_int | DNA_A_; }
      break;// 	C or G or T 	V
    case 'n': case 'N': case 'x': case 'X': case '?': {
      _fastSetNucleotideUnknown(c_pos);
      ambiguities.push_back((ambiguity_nucleotide_at_position)
      {{DNA_N_, 0.25, 0.25, 0.25, 0.25},c_pos});
      c_int = c_int | DNA_A_; }
      break;//  A or C or G or T  
    default:
      //if ( *c_str != '\0' )USER_WARNING("Incorrect input DNA string: \'" << *c_str << "\'");
      goto END_FOR_LOOP;
    } 
    //put the nucleotide into the int
    //c_int = c_int | n;


    //UPDATE
    //current position etc we are working on.
    const int next_c_pos = c_pos+1;

    // PENDING more  than 60 unneccessary conditions per b128
    //-----------------------
    //JUMPING INT
    if ( next_c_pos == next_int_jump_pos ){
      _mm_prefetch(c_str,_MM_HINT_T0);  
      c_b128 = set_int_b128(c_b128,c_int,INT_D_I(c_pos));
      
      //jump positions, i.e. charpositions that change b128 or int
      next_int_jump_pos = first_charpos_in_next_int(next_c_pos);
      
      c_int = 0;//same as get_int_b128(c_b128, INT_D_I(next_c_pos));

      //JUMP B128
      if ( next_c_pos == next_b128_jump_pos  ){
	set_b128(c_b128_ptr,c_b128);
	
	//jump to next b128 in the memory
	c_b128_ptr++;
	c_b128 = set_zero_b128();
	//jump positions, i.e. charpositions that change b128 or int
	next_b128_jump_pos = first_charpos_in_next_b128(next_c_pos);
      }
    }
    //-------------------------
    //update to the next position
    c_pos = next_c_pos;
    c_int = c_int << 2;
    
  }//end for
 END_FOR_LOOP:


  
  //LAST CHANGES INTO MEMORY
  //first justify all position in the int to the most significant side
  //print_bits_32(c_int);cout << endl;
  c_int = c_int << BIT_INT_I(c_pos);
  //print_bits_32(c_int);cout << endl;
  //  *c_b128_ptr = set_int_b128(c_b128,c_int,INT_D_I(c_pos));
  set_b128(c_b128_ptr,set_int_b128(c_b128,c_int,INT_D_I(c_pos)));
  

  //update global variables
  num_unknowns_ = _num_unknowns_;
  num_As_ = _num_As_;
  num_Cs_ = _num_Cs_;
  num_Gs_ = _num_Gs_;
  num_Ts_ = _num_Ts_;

  c_int = numChars-c_pos;
  numChars = c_pos;
  return c_int;
}




//-----------------------------------
// ACCESSING SPECIFIC NUCLEOTIDES
nucleotide
DNA_b128_String::getNucleotide(int pos) const{

  nucleotide n =  getNucleotideNOAMBIGUITY(pos);

  //check that not deleted or ambigious
  if ( n == DNA_UNKNOWN_ ){
    //check if ambigious (PENDING could do a binary search!)
    vector<ambiguity_nucleotide_at_position>::const_iterator iter = ambiguities.begin();
    for ( ; iter != ambiguities.end() && (*iter).position < pos ; ++iter ){
      //skip until pos
    }
    if ( iter != ambiguities.end() && (*iter).position == pos )
      return (*iter).ambiguity.n;
  }
  return n;
}
  
nucleotide
DNA_b128_String::getNucleotideNOAMBIGUITY(int pos) const{

  // _printPOS(pos);
  b128 c_b128 = get_b128(data + D_I(pos));
  
  int c_int = get_int_b128(c_b128, INT_D_I(pos));

  nucleotide n =  (nucleotide)((c_int >> BIT_INT_I(pos) ) & NUCLEOTIDE_INT_MASK);

  //check that not deleted or ambigious
  if ( n == DNA_A_ ){
    c_b128 = get_b128(unknownData + D_I(pos));
    c_int = get_int_b128(c_b128, INT_D_I(pos));
    if ( 0 != ((c_int >> BIT_INT_I(pos) ) & NUCLEOTIDE_INT_MASK) ){
      return DNA_UNKNOWN_;
    }
    //if not unknown or ambigious then it is an A
  }
  
  return n;
}

//returns the old nucleotide at the position
nucleotide
DNA_b128_String::setNucleotide(int pos, nucleotide n){

  b128 c_b128;
  int c_int;

  if ( n == DNA_NOT_ALLOWED ){
    cerr << "Can not set not allowed nucleotides" << endl;
    return DNA_NOT_ALLOWED;
  }
  
  //-----
  //update frequences
  nucleotide old = getNucleotide(pos);//PENDING slow
  switch ( old ){
  case DNA_A_:num_As_--;break;
  case DNA_C_:num_Cs_--;break;
  case DNA_G_:num_Gs_--;break;
  case DNA_T_:num_Ts_--;break;
  case DNA_UNKNOWN_:
    num_unknowns_--;
    //zero set unknown memory
    _fastClearNucleotideUnknown(D_I(pos),INT_D_I(pos),BIT_INT_I(pos));
    break;
  default://is ambigious remove it
    assert ( is_ambiguity_nucleotide(old) != 0 );
    _fastClearNucleotideUnknown(D_I(pos),INT_D_I(pos),BIT_INT_I(pos));
    //PENDING very slowDNA_b128_String!  (PENDING could do a binary search!)
    vector<ambiguity_nucleotide_at_position>::iterator iter = ambiguities.begin();
    for ( ; iter != ambiguities.end() ; ++iter ){
      if ( (*iter).position == pos ){
        ambiguities.erase(iter);        
        break;
      }
    }
  }

  //---
  //check the new nucleotide
  switch ( n ){
  case DNA_UNKNOWN_:
    num_unknowns_++;
    n = DNA_A_;
    //set unknown memory
    _fastSetNucleotideUnknown(D_I(pos),INT_D_I(pos),BIT_INT_I(pos));
    break;
  case DNA_A_:num_As_++;break;
  case DNA_C_:num_Cs_++;break;
  case DNA_G_:num_Gs_++;break;
  case DNA_T_:num_Ts_++;break;
  default:
    assert ( is_ambiguity_nucleotide(n) != 0 );
    _fastSetNucleotideUnknown(D_I(pos),INT_D_I(pos),BIT_INT_I(pos));

    vector<ambiguity_nucleotide_at_position>::iterator iter = ambiguities.begin();
    while ( iter != ambiguities.end() &&
            (*iter).position < pos )
      ++iter;
    ambiguity_nucleotide_at_position nap = {{n,0,0,0,0}, pos};
    ambiguities.insert(iter, nap );
    n = DNA_A_;
  }

  //----
  //set the nucleotide memory
  c_b128 = get_b128(data + D_I(pos));//PENDING SLOWDNA_b128_String! can save the D_I etc 
  c_int = get_int_b128(c_b128, INT_D_I(pos));
  //put 00 in the specific position
  c_int = c_int & (~(NUCLEOTIDE_INT_MASK << BIT_INT_I(pos)));
  //put 'n' at the specific position
  c_int = c_int | (n << BIT_INT_I(pos));

  //update the memory
  data[D_I(pos)] = set_int_b128(c_b128, c_int, INT_D_I(pos));

  return old;
}
  
// PRINTING
std::ostream &
DNA_b128_String::printOn(std::ostream &out) const{

  //PENDING VERY SLOW especially when there are ambiguities.
  for (size_t i = 0 ; i < numChars ; i++ )
    out << nucleotide2char(getNucleotide(i)) ;

  return out;
}

//--------------------------------------------------------
// EQUALS
bool operator== (const DNA_b128_String::ambiguity_nucleotide_at_position &a,  
		 const DNA_b128_String::ambiguity_nucleotide_at_position &b){
  
  return ( a.position != b.position ) ||
    (a.ambiguity.n != b.ambiguity.n);
}

bool 
DNA_b128_String::equals(const Object *o) const{
  const DNA_b128_String *s2 = (const DNA_b128_String *) o;
  if ( getNumChars() != s2->getNumChars() )
    return false;

  //PENDING can probably be improved by using the sse2 instructions for comparing
  if ( memcmp(data,s2->data,sizeof(b128)*getNumUsedDatas()) != 0 )
    return false;
  if ( memcmp(unknownData,s2->unknownData,sizeof(b128)*getNumUsedDatas()) != 0 )
    return false;

  return equal(ambiguities.begin(), ambiguities.end(), s2->ambiguities.begin());
}


//-------------------- AMBIGUITIES --------------------------------------------------------------

// RESOLVE AMBIGUITIES
void
DNA_b128_String::resolveAmbiguities(const DNA_b128_String &temp_str){

  std::vector<ambiguity_nucleotide_at_position>::iterator iter = ambiguities.begin();
  //  int numres = 0;
  for ( ; iter != ambiguities.end() ; ++iter ){
    ambiguity_nucleotide_at_position ambpos = *iter;
    nucleotide template_nucleotide = temp_str.getNucleotideNOAMBIGUITY(ambpos.position);

    //only resolve if the template is a regular nucleotide and if the template nucleotide is
    //contained in the ambiguity.    
    if ( template_nucleotide != DNA_UNKNOWN_ ){
      ambiguity_nucleotide tmp_amb = regularnucleotide2ambiguity_nucleotide(template_nucleotide);
      if ( is_ambiguity_contained(tmp_amb,ambpos.ambiguity) ){
	(*iter).ambiguity = tmp_amb;
      }
    }
  }
  //  cout << "total num amb " << getNumAmbiguities() << "    num resolv " << numres << endl;

}

void
DNA_b128_String::resolveAmbiguitiesUsingTransitionProbabilities(const DNA_b128_String &temp_str, 
								ML_string_distance mldist){

  std::vector<ambiguity_nucleotide_at_position>::iterator iter = ambiguities.begin();
 
  for ( ; iter != ambiguities.end() ; ++iter ){
    ambiguity_nucleotide_at_position ambpos = *iter;
    nucleotide template_nucleotide = temp_str.getNucleotideNOAMBIGUITY(ambpos.position);

    //only resolve if the template is a regular nucleotide and if the template nucleotide is
    //contained in the ambiguity.    
    if ( template_nucleotide != DNA_UNKNOWN_ ){
      ambiguity_nucleotide tmp_amb = regularnucleotide2ambiguity_nucleotide(template_nucleotide);
      if ( is_ambiguity_contained(tmp_amb,ambpos.ambiguity) ){
        ambiguity_nucleotide curramb = ambpos.ambiguity;
        float totalprob;
        switch( template_nucleotide ){
        case DNA_A_:
           totalprob = curramb.probA*mldist.A_A + curramb.probC*mldist.A_C + curramb.probG*mldist.A_G + curramb.probT*mldist.A_T;
          curramb.probA = curramb.probA*mldist.A_A/totalprob;
          curramb.probC = curramb.probC*mldist.A_C/totalprob;
          curramb.probG = curramb.probG*mldist.A_G/totalprob;
          curramb.probT = curramb.probT*mldist.A_T/totalprob;
          break;
        case DNA_C_:
           totalprob = curramb.probA*mldist.A_C + curramb.probC*mldist.C_C + curramb.probG*mldist.C_G + curramb.probT*mldist.C_T;
          curramb.probA = curramb.probA*mldist.A_C/totalprob;
          curramb.probC = curramb.probC*mldist.C_C/totalprob;
          curramb.probG = curramb.probG*mldist.C_G/totalprob;
          curramb.probT = curramb.probT*mldist.C_T/totalprob;
          break;
        case DNA_G_:
           totalprob = curramb.probA*mldist.A_G + curramb.probC*mldist.C_G + curramb.probG*mldist.G_G + curramb.probT*mldist.G_T;
          curramb.probA = curramb.probA*mldist.A_G/totalprob;
          curramb.probC = curramb.probC*mldist.C_G/totalprob;
          curramb.probG = curramb.probG*mldist.G_G/totalprob;
          curramb.probT = curramb.probT*mldist.G_T/totalprob;
          break;
        case DNA_T_:
           totalprob = curramb.probA*mldist.A_T + curramb.probC*mldist.C_T + curramb.probG*mldist.G_T + curramb.probT*mldist.T_T;
          curramb.probA = curramb.probA*mldist.A_T/totalprob;
          curramb.probC = curramb.probC*mldist.C_T/totalprob;
          curramb.probG = curramb.probG*mldist.G_T/totalprob;
          curramb.probT = curramb.probT*mldist.T_T/totalprob;
          break;
        default:
          PROG_ERROR("ERROR SHOULDNT COME HERE");
        }
        (*iter).ambiguity = curramb;
      }
    }
  }
}


//
// CORRECT DISTANCE COMPUTATION WITH AMBIGUITIES
//
simple_string_distance
DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(simple_string_distance sp,
                                                                         const DNA_b128_String &s1,
                                                                         const DNA_b128_String &s2){
  simple_string_distance real_distance = sp;

  //traverse the ambiguities and compute the correct distance.
  std::vector<ambiguity_nucleotide_at_position>::const_iterator i1 = s1.ambiguities.begin();
  std::vector<ambiguity_nucleotide_at_position>::const_iterator i2 = s2.ambiguities.begin();
  int pos1 = ( i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
  int pos2 = ( i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);

  //TRAVERSE
  // all ambiguities in order
  while( pos1 < INT_MAX || pos2 < INT_MAX ){

    //if two ambiguities match then there has been one then this has been
    //treated as a match of two unkowns.
    ambiguity_distance amdist;
    if ( pos1 == pos2 ){
      ambiguity_nucleotide an1 = (*i1).ambiguity;
      ambiguity_nucleotide an2 = (*i2).ambiguity;

      amdist = compute_ambiguity_distance(an1,an2);
      
      //---
      //update the distance
      real_distance.transitions += amdist.purine_transition_prob+amdist.pyrimidine_transition_prob;
      real_distance.transversions += amdist.transversion_prob;
      real_distance.deletedPositions -= 1;
      
      ++i1;
      ++i2;
      pos1 = ( i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
      pos2 = ( i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);
    }
    else if ( pos1 < pos2 ){
      ambiguity_nucleotide an1 = (*i1).ambiguity;
      nucleotide n = s2.getNucleotideNOAMBIGUITY(pos1);
      if ( n != DNA_UNKNOWN_ ){
        ambiguity_nucleotide an2 = regularnucleotide2ambiguity_nucleotide(n);
        amdist = compute_ambiguity_distance(an1,an2);
        
        //---
        //update the distance
        real_distance.transitions += amdist.purine_transition_prob+amdist.pyrimidine_transition_prob;
        real_distance.transversions += amdist.transversion_prob;
        real_distance.deletedPositions -= 1;
      }
      ++i1;
      pos1 = ( i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
    }
    else {// pos1 > pos2
      ambiguity_nucleotide an2 = (*i2).ambiguity;
      nucleotide n = s1.getNucleotideNOAMBIGUITY(pos2);
      if ( n != DNA_UNKNOWN_ ){
        ambiguity_nucleotide an1 = regularnucleotide2ambiguity_nucleotide(n);
        amdist = compute_ambiguity_distance(an1,an2);
        //---
        //update the distance
        real_distance.transitions += amdist.purine_transition_prob+amdist.pyrimidine_transition_prob;
        real_distance.transversions += amdist.transversion_prob;
        real_distance.deletedPositions -= 1;
      }
      
      ++i2;
      pos2 = ( i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);
    }

  }
  //---------------------
  
  return real_distance;  
}


TN_string_distance
DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(TN_string_distance sp,
                                                                         const DNA_b128_String &s1,
                                                                         const DNA_b128_String &s2){
  TN_string_distance real_distance = sp;

  //traverse the ambiguities and compute the correct distance.
  std::vector<ambiguity_nucleotide_at_position>::const_iterator i1 = s1.ambiguities.begin();
  std::vector<ambiguity_nucleotide_at_position>::const_iterator i2 = s2.ambiguities.begin();
  int pos1 = ( i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
  int pos2 = ( i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);

  //TRAVERSE
  // all ambiguities in order
  while( pos1 < INT_MAX || pos2 < INT_MAX ){

    //if two ambiguities match then there has been one then this has been
    //treated as a match of two unkowns.
    ambiguity_distance amdist;
    if ( pos1 == pos2 ){
      ambiguity_nucleotide an1 = (*i1).ambiguity;
      ambiguity_nucleotide an2 = (*i2).ambiguity;

      amdist = compute_ambiguity_distance(an1,an2);
      
      //---
      //update the distance
      real_distance.purine_transitions += amdist.purine_transition_prob;
      real_distance.pyrimidine_transitions += amdist.pyrimidine_transition_prob;
      real_distance.transversions += amdist.transversion_prob;
      real_distance.deletedPositions -= 1;
      
      ++i1;
      ++i2;
      pos1 = ( i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
      pos2 = ( i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);
    }
    else if ( pos1 < pos2 ){
      ambiguity_nucleotide an1 = (*i1).ambiguity;
      nucleotide n = s2.getNucleotideNOAMBIGUITY(pos1);
      if ( n != DNA_UNKNOWN_ ){
        ambiguity_nucleotide an2 = regularnucleotide2ambiguity_nucleotide(n);
        amdist = compute_ambiguity_distance(an1,an2);
        
        //---
        //update the distance
        real_distance.purine_transitions += amdist.purine_transition_prob;
        real_distance.pyrimidine_transitions += amdist.pyrimidine_transition_prob;
        real_distance.transversions += amdist.transversion_prob;
        real_distance.deletedPositions -= 1;
      }
      ++i1;
      pos1 = ( i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
    }
    else {// pos1 > pos2
      ambiguity_nucleotide an2 = (*i2).ambiguity;
      nucleotide n = s1.getNucleotideNOAMBIGUITY(pos2);
      if ( n != DNA_UNKNOWN_ ){
        ambiguity_nucleotide an1 = regularnucleotide2ambiguity_nucleotide(n);
        amdist = compute_ambiguity_distance(an1,an2);
        //---
        //update the distance
        real_distance.purine_transitions += amdist.purine_transition_prob;
        real_distance.pyrimidine_transitions += amdist.pyrimidine_transition_prob;
        real_distance.transversions += amdist.transversion_prob;
        real_distance.deletedPositions -= 1;
      }
      
      ++i2;
      pos2 = ( i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);
    }

  }
  //---------------------
  
  return real_distance;  
} 


simple_string_distance
DNA_b128_String::correctDistanceWithAmbiguitiesUsingTransitionProbabilities(simple_string_distance sp,
                                                                            ML_string_distance tp,
                                                                            const DNA_b128_String &s1,
                                                                            const DNA_b128_String &s2){
  
  simple_string_distance real_distance = sp;

  //traverse the ambiguities and compute the correct distance.
  std::vector<ambiguity_nucleotide_at_position>::const_iterator i1 = s1.ambiguities.begin();
  std::vector<ambiguity_nucleotide_at_position>::const_iterator i2 = s2.ambiguities.begin();
  int pos1 = ( i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
  int pos2 = ( i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);

  //TRAVERSE
  // all ambiguities in order
  while( pos1 < INT_MAX || pos2 < INT_MAX ){
    ambiguity_distance amdist;
    //if two ambiguities match then there has been one then this has been
    //treated as a match of two unkowns.
    if ( pos1 == pos2 ){
      ambiguity_nucleotide an1 = (*i1).ambiguity;
      ambiguity_nucleotide an2 = (*i2).ambiguity;

      amdist = compute_ambiguity_distance_using_transition_probabilities(an1,an2,tp);

      //---
      //update the distance
      real_distance.transitions += amdist.purine_transition_prob+amdist.pyrimidine_transition_prob;
      real_distance.transversions += amdist.transversion_prob;
      real_distance.deletedPositions -= 1;

      ++i1;
      ++i2;
      pos1 = ( i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
      pos2 = ( i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);
    }
    else if ( pos1 < pos2 ){
      ambiguity_nucleotide an1 = (*i1).ambiguity;
      nucleotide n = s2.getNucleotideNOAMBIGUITY(pos1);
      if ( n != DNA_UNKNOWN_ ){
	//        ambiguity_nucleotide an2 = regularnucleotide2ambiguity_nucleotide(n);
      
        amdist = compute_ambiguity_distance_using_transition_probabilities(n,an1,tp);
        //---
        //update the distance
        real_distance.transitions += amdist.purine_transition_prob+amdist.pyrimidine_transition_prob;
        real_distance.transversions += amdist.transversion_prob;
        real_distance.deletedPositions -= 1;
      }
      ++i1;
      pos1 = ( i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
    }
    else {// pos1 > pos2
      ambiguity_nucleotide an2 = (*i2).ambiguity;
      nucleotide n = s1.getNucleotideNOAMBIGUITY(pos2);
      if ( n != DNA_UNKNOWN_ ){
	//        ambiguity_nucleotide an1 = regularnucleotide2ambiguity_nucleotide(n);
      
        amdist = compute_ambiguity_distance_using_transition_probabilities(n,an2,tp);
        //---
        //update the distance
        real_distance.transitions += amdist.purine_transition_prob+amdist.pyrimidine_transition_prob;
        real_distance.transversions += amdist.transversion_prob;
        real_distance.deletedPositions -= 1;

      }
      
      ++i2;
      pos2 = ( i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);
    }

  }

  
  return real_distance;  
  
}


TN_string_distance
DNA_b128_String::correctDistanceWithAmbiguitiesUsingTransitionProbabilities(TN_string_distance sp,
                                                                            ML_string_distance tp,
                                                                            const DNA_b128_String &s1,
                                                                            const DNA_b128_String &s2){
  
  TN_string_distance real_distance = sp;

  //traverse the ambiguities and compute the correct distance.
  std::vector<ambiguity_nucleotide_at_position>::const_iterator i1 = s1.ambiguities.begin();
  std::vector<ambiguity_nucleotide_at_position>::const_iterator i2 = s2.ambiguities.begin();
  int pos1 = ( i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
  int pos2 = ( i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);

  //TRAVERSE
  // all ambiguities in order
  while( pos1 < INT_MAX || pos2 < INT_MAX ){
    ambiguity_distance amdist;
    //if two ambiguities match then there has been one then this has been
    //treated as a match of two unkowns.
    if ( pos1 == pos2 ){
      ambiguity_nucleotide an1 = (*i1).ambiguity;
      ambiguity_nucleotide an2 = (*i2).ambiguity;

      amdist = compute_ambiguity_distance_using_transition_probabilities(an1,an2,tp);
      //---
      //update the distance
      real_distance.purine_transitions += amdist.purine_transition_prob;
      real_distance.pyrimidine_transitions+= amdist.pyrimidine_transition_prob;
      real_distance.transversions += amdist.transversion_prob;
      real_distance.deletedPositions -= 1;

      ++i1;
      ++i2;
      pos1 = ( i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
      pos2 = ( i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);
    }
    else if ( pos1 < pos2 ){
      ambiguity_nucleotide an1 = (*i1).ambiguity;
      nucleotide n = s2.getNucleotideNOAMBIGUITY(pos1);
      if ( n != DNA_UNKNOWN_ ){
        ambiguity_nucleotide an2 = regularnucleotide2ambiguity_nucleotide(n);
      
        amdist = compute_ambiguity_distance_using_transition_probabilities(an1,an2,tp);
        //---
        //update the distance
        real_distance.purine_transitions += amdist.purine_transition_prob;
        real_distance.pyrimidine_transitions+= amdist.pyrimidine_transition_prob;
        real_distance.transversions += amdist.transversion_prob;
        real_distance.deletedPositions -= 1;
      }
      ++i1;
      pos1 = ( i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
    }
    else {// pos1 > pos2
      ambiguity_nucleotide an2 = (*i2).ambiguity;
      nucleotide n = s1.getNucleotideNOAMBIGUITY(pos2);
      if ( n != DNA_UNKNOWN_ ){
        ambiguity_nucleotide an1 = regularnucleotide2ambiguity_nucleotide(n);
      
        amdist = compute_ambiguity_distance_using_transition_probabilities(an1,an2,tp);
        //---
        //update the distance
        real_distance.purine_transitions += amdist.purine_transition_prob;
        real_distance.pyrimidine_transitions+= amdist.pyrimidine_transition_prob;
        real_distance.transversions += amdist.transversion_prob;
        real_distance.deletedPositions -= 1;

      }
      
      ++i2;
      pos2 = ( i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);
    }

  }
  //---------------------
  
  return real_distance;  
  
}







