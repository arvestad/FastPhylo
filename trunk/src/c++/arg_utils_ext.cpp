//--------------------------------------------------
//                                        
// File: argsutil_ext.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: arg_utils_ext.cpp,v 1.3 2006/08/23 07:40:08 isaac Exp $                                 
//
//--------------------------------------------------

#include "arg_utils_ext.hpp"
#include <cstdlib>

bool
get_boolean_option_value(int argc, char **argv, char *option_id){
  int i;

  i = has_option(argc,argv,option_id);
  if ( i ){
      return true;
  }
  return false;
}



bool
get_list_of_args(int argc, char **argv, char *option_id, std::vector<char*> &vec){

  int i = HAS_OPTION(option_id);
  if ( i <= 0 ) return false;
  i++;
  //read all file names that come after option_id until next option.
  for ( ; i < argc ; i++ ){
    if ( argv[i][0] != '-' ){ 
      vec.push_back(argv[i]);
    }
    else break;
  }

  return true;
}



bool
get_list_of_ints(int argc, char **argv, char *option_id, std::vector<int> &vec){

  int i = HAS_OPTION(option_id);
  if ( i <= 0 ) return false;
  i++;
  //read all file names that come after option_id until next option.
  for ( ; i < argc ; i++ ){
    if ( argv[i][0] != '-' ){ 
      vec.push_back(atoi(argv[i]));
    }
    else break;
  }

  return true;
}


bool
get_list_of_floats(int argc, char **argv, char *option_id, std::vector<float> &vec){

  int i = HAS_OPTION(option_id);
  if ( i <= 0 ) return false;
  i++;
  //read all file names that come after option_id until next option.
  for ( ; i < argc ; i++ ){
    if ( argv[i][0] != '-' ){ 
      vec.push_back(atof(argv[i]));
    }
    else break;
  }

  return true;
}
