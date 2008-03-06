
//--------------------------------------------------
//                                        
// File: argsutil.c                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: arg_utils.c,v 1.3 2006/12/19 12:20:52 isaac Exp $                                 
//
//--------------------------------------------------

#include "arg_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*  GET ARG  */
/*  goes through the arg array to find the id if the string */
/*  is found at index i then i is returned. if not found 0 is returned */

int
get_arg(int argc, char **argv, char *arg_id){
  int i;

  for ( i = 1 ; i < argc ; i++ ){
    if ( strcmp(argv[i],arg_id) == 0 )
      return i;
  }

  return 0;
}

/*  HAS OPTION */
/*  goes through the arg array searching for the option_id occuring as an option */
/*  returns the index of the occurrence or 0 if does not occurr. */

int
has_option(int argc, char **argv, char *option_id){
  
  int i;

  for ( i = 1 ; i < argc ; i++ ){
    if ( argv[i][0]=='-' ){
      if ( strcmp(argv[i],option_id) == 0 )
        return i;
    }
  }

  return 0;
}

/*  GET OPTION VALUE */
/*  searches for the option id and returns the value of the next variable */
/*  returns NULL if no value exists */

char *
get_option_value(int argc, char **argv, char *option_id){
  int i;

  i = has_option(argc,argv,option_id);
  if ( i ){
    if ( i + 1 < argc )
      return argv[i+1];
  }
  return NULL;
}

char *
set_int_option_value(int argc, char **argv, char *option_id, int *optvariable){

  char *chval = GET_OPTION_VAL(option_id);
  if ( chval != NULL )
    *optvariable = atoi(chval);

  return chval;
}

char *
set_float_option_value(int argc, char **argv, char *option_id, float *optvariable){

  char *chval = GET_OPTION_VAL(option_id);
  if ( chval != NULL )
    *optvariable = atof(chval);

  return chval;
}

char *
set_char_option_value(int argc, char **argv, char *option_id, char **optvariable){

  char *chval = GET_OPTION_VAL(option_id);
  if ( chval != NULL )
    *optvariable = chval;

  return chval;
}
  

/* /\*  SETUP STDIN/STDOUT *\/ */
/* /\*  Goes through searching for "stdin"/"stdout" and reopens the file  *\/ */
/* /\*  streams to the file names *\/ */

/* int */
/* setup_stdin(int argc, char **argv){ */
  
/*   char *file_name; */

/*   file_name = get_option_value(argc, argv, "stdin"); */
/*   if ( file_name != NULL ){ */
/*     stdin = freopen(file_name, "r",stdin); */
/*     return 1; */
/*   } */
/*   return 0; */
/* } */

/* int */
/* setup_stdout(int argc, char **argv){ */
  
/*   char *file_name; */

/*   file_name = get_option_value(argc, argv, "stdout"); */
/*   if ( file_name != NULL ){ */
/*     stdout = freopen(file_name, "w",stdout); */
/*     return 1; */
/*   } */
/*   return 0; */
/* } */

/* int */
/* cleanup_stdin(int argc, char **argv){ */
/*   char *file_name; */
/*   file_name = get_option_value(argc, argv, "stdin"); */
/*   if ( file_name != NULL ){ */
/*     fclose(stdin); */
/*     return 1; */
/*   } */
/*   return 0; */
/* } */


/* int */
/* cleanup_stdout(int argc, char **argv){ */
/*   char *file_name;  */
/*   file_name = get_option_value(argc, argv, "stdout"); */
/*   if ( file_name != NULL ){ */
/*     fclose(stdout); */
/*     return 1; */
/*   } */
/*   return 0; */
/* } */


/* int */
/* create_array_from_number_list(char *number_list, int **array){ */
/*   int num_numbers; */
/*   char *iter; */
/*   char list[1000]; */
/*   int i; */

/*   strcpy(list, number_list); */

/*   //count the number of ints */
/*   iter = list; */
/*   num_numbers = 1; */
/*   while ( *iter != '\0' ){ */
/*     if ( *iter == ',' ){ */
/*       num_numbers++; */
/*     } */
/*     iter++; */
/*   } */
  

/*   iter = strtok(list,","); */
/*   if ( iter == NULL ) */
/*     return -1; */
  
/*   //allocate the int array. */
/*   *array = (int *) malloc(num_numbers * sizeof(int)); */
/*   if ( *array == NULL ){ */
/*     fprintf(stderr,"OUT OF MEMORY\nfile: %s  func: %s  line: %d\n",__FILE__, __FUNCTION__, __LINE__); */
/*     exit(1); */
/*   } */
/*   (*array)[0] = atoi(iter); */
/*   i = 1; */
/*   while ( (iter = strtok(NULL,",")) != NULL ) { */
/*     (*array)[i] = atoi(iter); */
/*     i++; */
/*   }  */

/*   return num_numbers; */
/* } */






