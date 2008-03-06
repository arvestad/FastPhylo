//--------------------------------------------------
//                                        
// File: ambiguity_nucleotide.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: ambiguity_nucleotide.cpp,v 1.11 2006/12/27 13:40:49 isaac Exp $
//
//--------------------------------------------------

#include "ambiguity_nucleotide.hpp"



std::ostream&
operator<<(std::ostream & os, const ambiguity_distance &dist){
  os << "----" << std::endl;
  os << "purin ts " << dist.purine_transition_prob << std::endl;
  os << "pyrimidine ts " << dist.pyrimidine_transition_prob << std::endl;
  os << "transversion " << dist.transversion_prob << std::endl;
  os << "----" << std::endl;
  return os;
}

bool
is_ambiguity_contained(const nucleotide &n_, const nucleotide &set_){

  ambiguity_nucleotide n = nucleotide2ambiguity_nucleotide(n_);
  ambiguity_nucleotide set = nucleotide2ambiguity_nucleotide(set_);
  
  if ( n.probA > 0 && set.probA == 0 )
    return false;
  if ( n.probC > 0 && set.probC == 0 )
    return false;
  if ( n.probG > 0 && set.probG == 0 )
    return false;
  if ( n.probT > 0 && set.probT == 0 )
    return false;
  return true;
}










