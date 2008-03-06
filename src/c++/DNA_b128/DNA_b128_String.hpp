//--------------------------------------------------
//                                        
// File: DNA_b128_String.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: DNA_b128_String.hpp,v 1.16 2006/12/28 08:50:04 isaac Exp $                                 
//
//--------------------------------------------------

#ifndef DNA_b128_STRING_HPP
#define DNA_b128_STRING_HPP

#include <string>
#include "sse2_wrapper.h"
#include "nucleotide.hpp"
#include "ambiguity_nucleotide.hpp"
#include <iostream>
#include <vector>
#include "dna_pairwise_sequence_likelihood.hpp"
#include "Object.hpp"

//-----------------------------------------------
// DNA_b128_String
//
// This class contains a DNA sequences and allows for very fast
// computation of distances between two DNA sequences. This is
// achieved by using sse2 registers and an algorithm using bit
// fiddling and divide and conquer.
//
// The characters can be any of the nucleotides defined in nucleotide.h.
//
// EXAMPLE USAGE OF INIT
// std::string str1 = "agctagct-agctm"
// DNA_b128_String dna1(str1.length(),str1);
// std::string str2 = "ccctcccc-agctm"
// DNA_b128_String dna2(str2.length());
// dna2.append(str2);
//
// EXAMPLE USAGE OF DISTANCE COMPUTATION
// To compute the number of transitions and transversions use:
// simple_string_distance sd = DNA_b128_String::computeDistance(dna1,dna2);
//
// To compute the number of purine transition, pyrimidine transition, and transversions use:
// TN_string_distance tn = DNA_b128_String::computeTAMURANEIDistance(dna1,dna2);
//
// EXAMPLE USAGE OF INCLUDING AMBIGUITIES
// sd = DNA_b128_String::correctDistanceWithAmbiguitiesUsingTransitionProbabilities(sd,ml_dist,dna1,dna2);
// where ml_dist is a ML_string_distance with transition probabilities describing the Markov model.
//
//

class DNA_b128_String : public Object
{
public:
  //------------------ CONSTRUCTORS ------------------------------------
  //Creates an empty string which fits atleast capacity chars.
  DNA_b128_String();
  DNA_b128_String(int capacity);
 
  //COPY
  DNA_b128_String(const DNA_b128_String &str_b128);
  void operator=(const DNA_b128_String &str_b128);
 
  // Construct and fill atmost capacity chars.
  // See append() below.
  DNA_b128_String(int capacity, const std::string &str);
  DNA_b128_String(int capacity, const char *c_str);
  
  virtual ~DNA_b128_String();

  // If capacity>getTotalCapacity() then new memory is allocated.
  // Otherwise the chars and unknown datas are cleared.
  void reInitiate(size_t capacity=0);

  //---------------------------------------
  // CAPACITY
  // The total number of chars that can be fit into the string.
 size_t getTotalCapacity() const{
    return (_num_usable_datas()*64);
  }
  int getNumChars() const{
    return numChars;
  }
  //The number of more chars that can be appended to
  //the string.
  int getAdditionalCapacity() const{
    return (getTotalCapacity() - getNumChars());
  }
  
  //-----------------------------------
  // PRINTING 
  // Overridden from Object.
  // WARNINNG very slow.
  virtual std::ostream& printOn(std::ostream &out=std::cout) const;
    
  //-----------------------------------------
  // APPEND
  // Reads allowed nucleotide chars from str/in
  // and stops when a non DNA/RNA/unkown char is reached.
  // Spaces are just skipped. 
  //
  // Returns the number of chars read. 
  
  int append(const char *c_str);
  int append(const std::string &str){
    return append(str.c_str());
  }



  //-----------------------------------
  //ACCESSING SPECIFIC NUCLEOTIDES
  //WARNING: These are very slow! 
  nucleotide getNucleotide(int pos) const;
  //Assumes that the nucleotide is not an ambiguity, and returns
  //DNA_A_ if it is an ambiguity.
  nucleotide getNucleotideNOAMBIGUITY(int pos) const;
  //returns the old nucleotide at the position
  nucleotide setNucleotide(int pos, nucleotide n);

  //--------------------
  // EQUALS
  // VERY SLOW IMPLEMENTATION! 
  virtual bool equals(const Object *o) const;

  //------------------------------------
  // DISTANCE COMPUTATION
  // 
  // Computes the number of different missmatches between the two
  // strings.  All ambiguities are treated as unknowns. To update
  // distance computation to include the ambiguities call the
  // associated methods below.
  static simple_string_distance computeDistance(const DNA_b128_String &s1,
                                                const DNA_b128_String &s2);
  static TN_string_distance computeTAMURANEIDistance(const DNA_b128_String &s1,
                                                     const DNA_b128_String &s2);


  //the number of gaps is the number of missmatches of '-' toward something that is not '-'.
  static int getNumberOfGaps(const DNA_b128_String &s1,
                             const DNA_b128_String &s2,
                             simple_string_distance d){
    return  d.deletedPositions - (s1.num_unknowns_ + s2.num_unknowns_ - d.deletedPositions);
  }


  //------------------------ AMBIGUITIES -------------------------------------------------------

  //HAS AMBIGUITY
  bool hasAmbiguities() const{
    return (ambiguities.size() >0);
  }
  int getNumAmbiguities() const{
    return ambiguities.size();
  }  

  // RESOLVE AMBIGUITIES
  // Uses the template string to resolve ambiguities. If the template string has an ambiguity or a gap
  // in one of the positions that the current string has then the position is remains unresolved.
  void resolveAmbiguities(const DNA_b128_String &temp_str);
  void resolveAmbiguitiesUsingTransitionProbabilities(const DNA_b128_String &temp_str, ML_string_distance mldist);

  // CORRECT DISTANCE WITH AMBIGUITIES
  //
  // There are two models for this:
  //      1) use the frequencies of nucleotides to assign probabilities of each ambiguity (uniform or background).
  //      2) use the frequencies AND the transprobabilities to compute the ambiguity distance
  
  //Before this function is called make sure that calcAmbiguityProbabilities() have been called.
  static simple_string_distance correctDistanceWithAmbiguitiesUsingBackgroundFrequences(simple_string_distance sp,
                                                                                        const DNA_b128_String &s1,
                                                                                        const DNA_b128_String &s2);
  static TN_string_distance correctDistanceWithAmbiguitiesUsingBackgroundFrequences(TN_string_distance sp,
                                                                                    const DNA_b128_String &s1,
                                                                                    const DNA_b128_String &s2);

  //Takes a compute distance in which the ambiguities have not been considered and transition probabilities.
  //And computes a new distance in which the ambiguities have been considered.
  static simple_string_distance correctDistanceWithAmbiguitiesUsingTransitionProbabilities(simple_string_distance sp,
                                                                                           ML_string_distance tp,
                                                                                           const DNA_b128_String &s1,
                                                                                           const DNA_b128_String &s2);
  static TN_string_distance correctDistanceWithAmbiguitiesUsingTransitionProbabilities(TN_string_distance sp,
                                                                                           ML_string_distance tp,
                                                                                           const DNA_b128_String &s1,
                                                                                           const DNA_b128_String &s2);

  //-----------------------------
  // SET AMBIGUITY DISTRIBUTION
  //
  //These functions compute the probabilities of every ambiguity symbol.
  //E.g. 'M' is either an 'A' or a 'C' and if basefreqA=3 and basefreqC=2
  //then the probability of seeing an 'A' is 3/5 and for 'C' 2/5.
  void calcAmbiguityProbabilities(int basefreqA, int basefreqC,
                                  int basefreqG, int basefreqT){
    std::vector<ambiguity_nucleotide_at_position>::iterator iter = ambiguities.begin();
    for ( ; iter != ambiguities.end() ; ++iter ){
      (*iter).ambiguity = nucleotide2ambiguity_nucleotide((*iter).ambiguity.n, basefreqA, basefreqC, basefreqG, basefreqT);
    }
  }
  //set the distribution of each ambiguity to be uniform over the allowed nucleotides in the position.
  void calcAmbiguityProbabilitiesUNIFORM(){
    std::vector<ambiguity_nucleotide_at_position>::iterator iter = ambiguities.begin();
    for ( ; iter != ambiguities.end() ; ++iter ){
      (*iter).ambiguity = nucleotide2ambiguity_nucleotideUNIFORM((*iter).ambiguity.n);
    }
  }  


  //---------------------------------------------------------------------
  // BASE FREQUENCES
  typedef struct{
    int num_As_;
    int num_Cs_;
    int num_Gs_;
    int num_Ts_;
    int num_unknowns_;
    int num_ambiguities_;
  } base_frequences;
  
  base_frequences getBaseFrequences() const{
    base_frequences bf ={num_As_,num_Cs_,num_Gs_,num_Ts_,num_unknowns_,ambiguities.size()};
    return bf;
  }
  
  //---------------------------- PRIVATE -----------------------------------

private:
  // THE DATA
  // The sequences is stored in three different types:
  //
  // 1) data: contains the nucleotides. Each nucleotide is made up of
  // two bits as described in nucleotide.h. If a position contains a
  // gap '-' or an ambiguity nucleotide then the two bits in data is
  // an DNA_A_, i.e. 00.  To simplify reading the chars in each b128
  // are in reversed order. For example, char at position 0 is located
  // at bits 126 and 127 of the first b128 and char at position 127 is
  // located at bits 0 and 1 of the second b128.
  //
  //
  // 2) unknownData: Is of the same size as data. If a position
  // contains a gap or an ambiguity nucleotide then the two associated
  // bits are 11, otherwise 00. A position contains a gap if the
  // ambiguities vector does not contain an ambiguity at the position.
  //
  // 3) ambiguities: is a vector containing all ambiguity nucleotides
  // and their positions. WARNING. Note that this implementation
  // depends on that there are very few ambiguities. Especially the
  // get/set methods are very slow due to this array.

  size_t numChars;    
  b128 *data;
  size_t numDatas;
  b128 *unknownData;
  
  typedef struct{
    ambiguity_nucleotide ambiguity;
    int position;
  } ambiguity_nucleotide_at_position;
  
  std::vector<ambiguity_nucleotide_at_position> ambiguities;
  friend bool operator== (const DNA_b128_String::ambiguity_nucleotide_at_position &a,  
			  const DNA_b128_String::ambiguity_nucleotide_at_position &b);
  
  size_t _num_usable_datas() const {
    return (numDatas -1);//the last b128 is not used!
  }
  //the number of chars currently in the string.

  int getNumUsedDatas() const {
    return ((numChars%64) == 0 ? numChars/64 : (numChars/64 +1));
  }
  
  //Base frequences (neccesary in the TN93 correction formula)
  size_t num_As_, num_Cs_, num_Gs_, num_Ts_;
  size_t num_unknowns_;

  //
  // Allocates and frees the memory. To change the capacity of a
  // sequence call void DNA_b128_String::reInitiate(size_t capacity).
  //
  void _init_mem(int capacity);
  void _free_mem();

  //
  // Sets and clears the position in unkownData
  //
  void _fastSetNucleotideUnknown(int data_pos,
                                 int int_in_data_pos,
                                 int char_in_int_pos){ 
    b128 c_b128 = get_b128(unknownData + data_pos);
    int c_int = get_int_b128(c_b128, int_in_data_pos);
    c_int = c_int | (NUCLEOTIDE_INT_MASK << char_in_int_pos);
    set_b128(unknownData+data_pos,set_int_b128(c_b128, c_int, int_in_data_pos));
  }
  void _fastSetNucleotideUnknown(int charpos){
    b128 *u = unknownData + D_I(charpos);
    b128 c_b128 = get_b128(u);
    const int intpos=INT_D_I(charpos);
    int c_int = get_int_b128(c_b128, intpos);
    c_int = c_int | (NUCLEOTIDE_INT_MASK << BIT_INT_I(charpos));
    set_b128(u,set_int_b128(c_b128, c_int, intpos));
  }

  void _fastClearNucleotideUnknown(int data_pos,
                                   int int_in_data_pos,
                                   int char_in_int_pos){ 
    b128 c_b128 = get_b128(unknownData + data_pos);
    int c_int = get_int_b128(c_b128, int_in_data_pos);
    c_int = c_int ^ (NUCLEOTIDE_INT_MASK << char_in_int_pos);//clears
    set_b128(unknownData+data_pos,set_int_b128(c_b128, c_int, int_in_data_pos));
  }
  
  //-----------------------------------
  //MEMORY ACCESS
  //the data index
  static int D_I(int charpos){
    return (charpos >> 6);// (charpos/64)
  }
  //the index of the int in the b128
  static int INT_D_I(int charpos){
    return (3 - ( (charpos & 0x3f)>> 4)); //3-((charpos%64)/16)
  }
  //the index of the bit in the int
  //(charpos=0 is located at bit index 30 and 31)
  static int BIT_INT_I(int charpos){
    return (30 - ((charpos & 0xf)<<1)); //30- 2*(charpos%16)
  }
  //if the charpos is the in a b128 slot
  static bool first_char_in_dataslot(int charpos){
    return ((charpos & 0x3f) == 0);
  }
  //if the charpos is the in a b128 slot
  static bool first_char_in_intslot(int charpos){
    return ((charpos & 0xf) == 0);
  }
  //returns the next charpos which is in another b128
  static int first_charpos_in_next_b128(int charpos){
    return (((charpos >> 6)+1)<<6);
  }//returns the next charpos which is in another b128
  static int first_charpos_in_next_int(int charpos){
    return  (((charpos >> 4)+1)<<4);
  }

  //utility
  static void _printPOS(int charpos){
    std::cout << "charpos: " << charpos << "    D_I: " << 
      D_I(charpos) << "  INT_D_I: " <<INT_D_I(charpos) << "   BIT_INT_I: " << BIT_INT_I(charpos) << std::endl; 
  }
  
};


#endif // DNA_b128_STRING_HPP




