//--------------------------------------------------
//                                        
// File: argsutil_ext.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: arg_utils_ext.hpp,v 1.4 2006/12/03 14:01:25 isaac Exp $                                 
//
//
// Provides extended functionality to command line argument handeling in C++
//--------------------------------------------------
#ifndef EXTARGSUTIL_HPP
#define EXTARGSUTIL_HPP

#include <vector>
#include "arg_utils.h"
#include <stdlib.h>


/*  GET BOOLEAN OPTION VALUE */
/*  searches for the option id and returns true if it exists otherwise false.*/
#define GET_BOOLEAN_OPTION_VAL(opt_id) get_boolean_option_value(argc,argv,opt_id)
bool
get_boolean_option_value(int argc, char **argv, char *option_id);


/* GET LIST OF ARGS */
/* Returns false if opt id is not found otherwise it returns true and fills in the vector. */
/* the vector is filled in with all values following the opt id until the next id*/
#define GET_LIST_OF_ARGS(opt_id,char_vec) get_list_of_args(argc,argv,opt_id,char_vec)
bool
get_list_of_args(int argc, char **argv, char *option_id, std::vector<char*> &vec);


#define GET_LIST_OF_INTS(opt_id,char_vec) get_list_of_ints(argc,argv,opt_id,char_vec)
bool
get_list_of_ints(int argc, char **argv, char *option_id, std::vector<int> &vec);

#define GET_LIST_OF_FLOATS(opt_id,char_vec) get_list_of_floats(argc,argv,opt_id,char_vec)
bool
get_list_of_floats(int argc, char **argv, char *option_id, std::vector<float> &vec);




#endif // EXTARGSUTIL_HPP
