//--------------------------------------------------
//                                        
// File: Sequences2DistanceMatrix.hpp
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: Sequences2DistanceMatrix.hpp,v 1.7 2006/12/25 18:25:50 isaac Exp $                                 
//
//--------------------------------------------------

#ifndef SEQUENCES2DISTANCEMATRIX_HPP
#define SEQUENCES2DISTANCEMATRIX_HPP

#include "DistanceMatrix.hpp"
#include <string>
#include "DNA_b128_String.hpp"
#include <fstream>
#include "dna_pairwise_sequence_likelihood.hpp"
#include "SequenceTree.hpp"

//
// This file contains functions for creating a distance matrix from
// DNA sequences using different models of evolution.
//  
// The general function is:
//
// void fillMatrix(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs,
//		  sequence_translation_model trans_model);
//
// It fills the distance matrix dm using the input sequences and the
// translation model. Note that the identifiers of the rows in the
// matrix are untouched and need to be set by the caller.
//

//This struct describes how the distance should be computed from the strings.
typedef struct {

  //the correction formula to use
  sequence_model model;
  //if ambiguities should be included
  bool no_ambiguities;
  //if nearests neighbor resulotion should be used
  bool no_ambig_resolve;
  //if the MC trnasition probabilities should  be used for the ambiguities
  bool no_transition_probs;
  //if the base frequences should be used for the ambiguities.
  bool use_base_freqs;
  //FIXED TS TV RATIO
  float tstvratio;
  float pyrtvratio;
  bool no_tstvratio;
} sequence_translation_model;



//---------------------------------------------------------------
//Fills the distance matrix according to the sequence_translation_model.

void fillMatrix(StrDblMatrix &dm, std::vector<DNA_b128_String> &b128_strings,
		sequence_translation_model trans_model);


//--------------------------------------------------------------------
//Functions for specific models.
void 
fillMatrix_Hamming(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs,
		   sequence_translation_model trans_model);

void fillMatrix_JC(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs,
		   sequence_translation_model trans_model);

void fillMatrix_K2P(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs,
		    sequence_translation_model trans_model);
  


void fillMatrix_TN93(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs, 
		     DNA_b128_String::base_frequences freqs,
		     sequence_translation_model trans_model);

//-----------




//---------------------------------------------
// INSTANTIATION OF DNA_B128_Strings
//
//

//Takes a vector of Sequences and converts them to b128 strings.
void
Sequences2DNA_b128(std::vector<Sequence> &seqs, std::vector<DNA_b128_String> &b128); 

//Reads a phylip file and instantiates a vector of b128 strings. The
//names of the sequences are stored in names s.t. names[i] is the name
//of b128_strings[i].
void DNA_b128_StringsFromPHYLIP(std::ifstream &fin, std::vector<std::string> &names, std::vector<DNA_b128_String> &b128_strings);

//
// Creates a bootstrapped set of b128_strings from the input sequences.
//
void 
bootstrapSequences(const std::vector<Sequence> &seqs, std::vector<DNA_b128_String> &b128_strings);
  
#endif // SEQUENCES2DISTANCEMATRIX_HPP











