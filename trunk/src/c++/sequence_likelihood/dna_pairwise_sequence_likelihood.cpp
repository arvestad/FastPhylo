//--------------------------------------------------
//                                        
// File: dna_pairwise_sequence_likelihood.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: dna_pairwise_sequence_likelihood.cpp,v 1.1.1.1 2006/01/23 23:58:45 isaac Exp $                                 
//
//--------------------------------------------------

#include "dna_pairwise_sequence_likelihood.hpp"


#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

simple_string_distance
convert_TN_string_distance_to_simple(TN_string_distance &d){

  simple_string_distance sd;

  sd.deletedPositions = d.deletedPositions;
  sd.transitions = d.purine_transitions + d.pyrimidine_transitions;
  sd.transversions = d.transversions;

  return sd;
}


bool
operator==(const simple_string_distance &d1, const simple_string_distance &d2){
  return 
    (d1.deletedPositions == d2.deletedPositions) &&
    (d1.transitions == d2.transitions) &&
    (d1.transversions == d2.transversions);
}
bool
operator==(const TN_string_distance &d1, const TN_string_distance &d2){
  return 
    (d1.deletedPositions == d2.deletedPositions) &&
    (d1.purine_transitions == d2.purine_transitions) &&
    (d1.pyrimidine_transitions == d2.pyrimidine_transitions) &&
    (d1.transversions == d2.transversions);
}

std::ostream&
operator<<(std::ostream & os, const simple_string_distance &d){
  os << "{ deletedPositions = " << d.deletedPositions << " , transitions = " << d.transitions <<" , transversions = " << d.transversions << " }";
  return os;
}


std::ostream&
operator<<(std::ostream & os, const TN_string_distance &d){
  os << "{ deletedPositions = " << d.deletedPositions << " , purine_trans = " << d.purine_transitions << " , pyrimidine_trans = " << d.pyrimidine_transitions <<" , transversions = " << d.transversions << " }";
  return os;
}

std::ostream&
operator<<(std::ostream & os, const ML_string_distance &td){
  os << "----" << endl;
  os << "Distance = " << td.distance << endl;
  os << "Probability of observing respective changes:" << endl;
  os << "  " <<  "    A     " << "    C     " << "    G     " << "    T    " << "\n";
  os << "A " << setw(10) << td.A_A << setw(10) << td.A_C << setw(10) << td.A_G << setw(10) << td.A_T << "\n";
  os << "C " << "          " << setw(10) << td.C_C << setw(10) << td.C_G << setw(10) << td.C_T << "\n";
  os << "G " << "          " << "          " << setw(10) << td.G_G << setw(10) << td.G_T << "\n";
  os << "T " << "          " << "          " << "          " << setw(10) << td.T_T << endl;
  os << "----" << endl;
  return os;
}
