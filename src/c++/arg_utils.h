
//--------------------------------------------------
//                                        
// File: argsutil.h                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: arg_utils.h,v 1.3 2006/12/19 12:20:52 isaac Exp $                                 
//
//--------------------------------------------------

#ifndef ARGSUTIL_H
#define ARGSUTIL_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

  
#define STREQ(s1,s2) (strcmp(s1,s2)==0)  
/*  GET ARG  */
/*  goes through the arg array to find the id if the string */
/*  is found at index i then i is returned. if not found 0 is returned */
#define GET_ARG(arg_id) get_arg(argc,argv,arg_id)
  
int
get_arg(int argc, char **argv, char *arg_id);

/*  HAS OPTION */
/*  goes through the arg array searching for the option_id occuring as an option */
/*  returns the index of the occurrence or 0 if does not occurr. */
#define HAS_OPTION(opt_id) has_option(argc,argv,opt_id)
int
has_option(int argc, char **argv, char *option_id);
  

/*  GET OPTION VALUE */
/*  searches for the option id and returns the value of the next variable */
/*  returns NULL if no value exists */
#define GET_OPTION_VAL(opt_id) get_option_value(argc,argv,opt_id)
char *
get_option_value(int argc, char **argv, char *option_id);
 
/* GETSET */
/* Checks for the option and sets the option variable if it exists. */
/* Returns the char pointer to the option value*/
#define SET_INT_OPTION_VAL(opt_id,opt_var) set_int_option_value(argc,argv,opt_id,opt_var)
char *
set_int_option_value(int argc, char **argv, char *option_id, int *optvariable);

#define SET_FLOAT_OPTION_VAL(opt_id,opt_var) set_float_option_value(argc,argv,opt_id,opt_var)
char *
set_float_option_value(int argc, char **argv, char *option_id, float *optvariable);

#define SET_CHAR_OPTION_VAL(opt_id,opt_var) set_char_option_value(argc,argv,opt_id,opt_var)
char *
set_char_option_value(int argc, char **argv, char *option_id, char **optvariable);

  
/* /\*  SETUP STDIN/STDOUT *\/ */
/* /\*  Goes through searching for "stdin"/"stdout" and reopens the file  *\/ */
/* /\*  streams to the file names *\/ */

/* int */
/* setup_stdin(int argc, char **argv); */

/* int */
/* setup_stdout(int argc, char **argv); */

/* /\*  CLEANUP STDIN/STDOUT *\/ */
/* /\*  closes stdin/stdout if it exists in the args *\/ */

/* int */
/* cleanup_stdin(int argc, char **argv); */

/* int */
/* cleanup_stdout(int argc, char **argv); */

/* /\* NUMBER LIST *\/ */
/* /\* takes comma seperated list of numbers "1,2,3" allocates*\/ */
/* /\* an array and returns the number of elements (-1 if error) *\/ */

/* int */
/* create_array_from_number_list(char *number_list, int **array); */

#ifdef __cplusplus
}
#endif  

  
#endif /*ARGSUTIL.H*/








