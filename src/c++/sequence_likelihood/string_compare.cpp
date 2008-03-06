//--------------------------------------------------
//                                        
// File: string_compare.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: string_compare.cpp,v 1.4 2006/12/08 11:09:14 isaac Exp $                                 
//
//--------------------------------------------------

#include "string_compare.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

#include "nucleotide.hpp"
#include "ambiguity_nucleotide.hpp"
#include "log_utils.hpp"

using namespace std;


int hamming_distance(const std::string &s1,
		     const std::string &s2){

  int ham = 0;
  if( s1.length()!=s2.length() ){
    USER_ERROR("Sequences of different length " << s1.length() << "!=" << s2.length() << ".");
  }
  for ( size_t i = 0 ; i < s1.length() ; i++ ){
    //cout << s1[i] << " ? " << s2[i] <<" , ";
    if ( s1[i] != s2[i] )
      ham++;
  }
  
  return ham;
}




TN_string_distance
TN_string_compare(const std::string &s1,
                  const std::string &s2){

  float PURTS = 0;
  float PYRTS = 0;
  float TV = 0;
  int deleted = 0;
  int matchDeletions = 0;
  
  for ( unsigned int i = 0 ; i < s1.length() ; i++ ){
    nucleotide n1 = char2nucleotide(s1[i]);
    nucleotide n2 = char2nucleotide(s2[i]);
    nucleotide_difference diff = get_nucleotide_difference(n1,n2);
    switch( diff ){
    case EQUAL: break;//do nothing
    case DELETION: deleted++; if ( n1==n2 ) matchDeletions++; break;
    case PURINE_TRANSITION: PURTS++; break;
    case PYRIMIDINE_TRANSITION: PYRTS++; break;
    case TRANSVERSION: TV++; break;
    case AMBIGIOUS:
      USER_ERROR("amb in naive not handled");
    }
  }

  TN_string_distance d = {  deleted, PURTS,PYRTS, TV};
  return d;
}



int
complete_dna_string_compare(float divergence_matrix[4][4],//a 4x4 matrix will be filled in with the frequences
                            const std::string &s1,
                            const std::string &s2){

  vector<int> ambig_pos;
  int deleted = 0;
  float tmp_matrix[4][4];
  
  for ( int i = 0 ; i < 4 ; i++)
    for ( int j = 0 ; j < 4 ; j++ ){
      divergence_matrix[i][j] = 0;
      tmp_matrix[i][j] = 0;
    }
  for ( unsigned int  i = 0 ; i < s1.length() ; i++ ){
    nucleotide n1 = char2nucleotide(s1[i]);
    nucleotide n2 = char2nucleotide(s2[i]);
    nucleotide_difference diff = get_nucleotide_difference(n1,n2);
    switch( diff ){
    case DELETION: deleted++; break;
    case AMBIGIOUS:
      ambig_pos.push_back(i);
      break;
    default:
      divergence_matrix[n1][n2]++;
    }
  }

  //  cout << "div matrix before update" << endl;
  // print_divergence_matrix(divergence_matrix);
  
  //update ambiguities according to Swofford
  vector<int>::iterator iter = ambig_pos.begin();
  for ( ; iter != ambig_pos.end() ; ++iter ){
    int pos = *iter;
    nucleotide n1 = char2nucleotide(s1[pos]);
    nucleotide n2 = char2nucleotide(s2[pos]);
    ambiguity_nucleotide a1 = nucleotide2ambiguity_nucleotide(n1);
    ambiguity_nucleotide a2 = nucleotide2ambiguity_nucleotide(n2);


    float total = 0;
    for ( int i = 0 ; i <= 3 ; i++ ){
        ambiguity_nucleotide ai = nucleotide2ambiguity_nucleotide((nucleotide)i);
      for ( int j = 0 ; j <= 3 ; j++ ){
        ambiguity_nucleotide aj = nucleotide2ambiguity_nucleotide((nucleotide)j);
        if ( is_ambiguity_contained(ai,a1) && is_ambiguity_contained(aj,a2) ){
          total += divergence_matrix[i][j] + divergence_matrix[j][i];//PENDING
        }
      }
    }
    for ( int i = 0 ; i <= 3 ; i++ ){
      for ( int j = 0 ; j <= 3 ; j++ ){
        ambiguity_nucleotide ai = nucleotide2ambiguity_nucleotide((nucleotide)i);
        ambiguity_nucleotide aj = nucleotide2ambiguity_nucleotide((nucleotide)j);
        if ( is_ambiguity_contained(ai,a1) && is_ambiguity_contained(aj,a2) && total != 0){
          tmp_matrix[i][j] += (divergence_matrix[i][j] + divergence_matrix[j][i])/total;//PENDING
        }
      }
    }

  }//end for pos

  for ( int i = 0 ; i <= 3 ; i++ )
    for ( int j = 0 ; j <=3 ; j++ )
      divergence_matrix[i][j] = (divergence_matrix[i][j]+tmp_matrix[i][j]);

  // cout << "div matrix after update" << endl;
  // print_divergence_matrix(divergence_matrix);
  
  return deleted;
}



TN_string_distance
divergence_matrix_2_TN_distance(float divergence_matrix[4][4], int deleted){

  float pu = divergence_matrix[DNA_A_][DNA_G_] +
    divergence_matrix[DNA_G_][DNA_A_];
  float py = divergence_matrix[DNA_C_][DNA_T_] +
    divergence_matrix[DNA_T_][DNA_C_];
  float tv = divergence_matrix[DNA_A_][DNA_C_] +
    divergence_matrix[DNA_C_][DNA_A_] +
    divergence_matrix[DNA_G_][DNA_C_] +
    divergence_matrix[DNA_C_][DNA_G_] +
    divergence_matrix[DNA_A_][DNA_T_] +
    divergence_matrix[DNA_T_][DNA_A_] +
    divergence_matrix[DNA_G_][DNA_T_] +
    divergence_matrix[DNA_T_][DNA_G_];
  
    
  TN_string_distance tn = {deleted,pu,py,tv};

  return tn;
}


void
print_divergence_matrix(std::ostream &out, float div_matrix[4][4]){


  out << "     ";
  for ( int i = 0 ; i < 4 ; i++ ){
    out << setw(10) << nucleotide2char((nucleotide)i);
  }
  out << endl;

  for ( int i = 0 ; i < 4 ; i++ ){
    out <<  nucleotide2char((nucleotide)i) << "    ";
    for ( int j = 0 ; j < 4 ; j++ ){
      out << setw(10) << div_matrix[i][j];
    }
    out << endl;
  }
}













