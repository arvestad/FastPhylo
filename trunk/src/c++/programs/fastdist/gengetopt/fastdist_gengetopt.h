/** @file fastdist_gengetopt.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef FASTDIST_GENGETOPT_H
#define FASTDIST_GENGETOPT_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name */
#define CMDLINE_PARSER_PACKAGE "fastdist"
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION "0.9"
#endif

enum enum_input_format { input_format_arg_phylip_multialignment = 0 , input_format_arg_xml };
enum enum_output_format { output_format_arg_phylip_dm = 0 , output_format_arg_xml };
enum enum_distance_function { distance_function_arg_JC = 0 , distance_function_arg_K2P, distance_function_arg_FAKE_F84, distance_function_arg_TN93, distance_function_arg_HAMMING };
enum enum_ambiguity_frequency_model { ambiguity_frequency_model_arg_UNI = 0 , ambiguity_frequency_model_arg_BASE };
enum enum_method { method_arg_NJ = 0 , method_arg_FNJ, method_arg_BIONJ };

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  char * outfile_arg;	/**< @brief output filename. If not specifed, output is written to stdout.  */
  char * outfile_orig;	/**< @brief output filename. If not specifed, output is written to stdout original value given at command line.  */
  const char *outfile_help; /**< @brief output filename. If not specifed, output is written to stdout help description.  */
  enum enum_input_format input_format_arg;	/**< @brief input format (default='phylip_multialignment').  */
  char * input_format_orig;	/**< @brief input format original value given at command line.  */
  const char *input_format_help; /**< @brief input format help description.  */
  enum enum_output_format output_format_arg;	/**< @brief output format (default='phylip_dm').  */
  char * output_format_orig;	/**< @brief output format original value given at command line.  */
  const char *output_format_help; /**< @brief output format help description.  */
  enum enum_distance_function distance_function_arg;	/**< @brief Distance function (default='K2P').  */
  char * distance_function_orig;	/**< @brief Distance function original value given at command line.  */
  const char *distance_function_help; /**< @brief Distance function help description.  */
  int bootstraps_arg;	/**< @brief Bootstrap num times and create matrix for each (default='0').  */
  char * bootstraps_orig;	/**< @brief Bootstrap num times and create matrix for each original value given at command line.  */
  const char *bootstraps_help; /**< @brief Bootstrap num times and create matrix for each help description.  */
  int no_incl_orig_flag;	/**< @brief If the distance matrix from the original sequences should be included (default=off).  */
  const char *no_incl_orig_help; /**< @brief If the distance matrix from the original sequences should be included help description.  */
  int seed_arg;	/**< @brief Random seed. If not specified the current timestamp will be used.  */
  char * seed_orig;	/**< @brief Random seed. If not specified the current timestamp will be used original value given at command line.  */
  const char *seed_help; /**< @brief Random seed. If not specified the current timestamp will be used help description.  */
  int no_ambiguities_flag;	/**< @brief Ignore ambiguities (default=off).  */
  const char *no_ambiguities_help; /**< @brief Ignore ambiguities help description.  */
  int no_ambig_resolve_flag;	/**< @brief Specifies that ambigious symbols should not be resolved by nearest neighbor (default=off).  */
  const char *no_ambig_resolve_help; /**< @brief Specifies that ambigious symbols should not be resolved by nearest neighbor help description.  */
  int no_transprob_flag;	/**< @brief Specifies that the transition probabilities should not be used in the ambiguity model (default=off).  */
  const char *no_transprob_help; /**< @brief Specifies that the transition probabilities should not be used in the ambiguity model help description.  */
  enum enum_ambiguity_frequency_model ambiguity_frequency_model_arg;	/**< @brief Ambiguity frequency model (default='UNI').  */
  char * ambiguity_frequency_model_orig;	/**< @brief Ambiguity frequency model original value given at command line.  */
  const char *ambiguity_frequency_model_help; /**< @brief Ambiguity frequency model help description.  */
  float tstvratio_arg;	/**< @brief Transition/transvertion ratio for purine transitions ( for the TN model ) (default='2.0').  */
  char * tstvratio_orig;	/**< @brief Transition/transvertion ratio for purine transitions ( for the TN model ) original value given at command line.  */
  const char *tstvratio_help; /**< @brief Transition/transvertion ratio for purine transitions ( for the TN model ) help description.  */
  float pyrtvratio_arg;	/**< @brief Transition/transvertion ratio for  pyrimidines transitions ( for the TN model ) (default='2.0').  */
  char * pyrtvratio_orig;	/**< @brief Transition/transvertion ratio for  pyrimidines transitions ( for the TN model ) original value given at command line.  */
  const char *pyrtvratio_help; /**< @brief Transition/transvertion ratio for  pyrimidines transitions ( for the TN model ) help description.  */
  int no_tstvratio_flag;	/**< @brief If given fixed ts/tv ratios will not be used (default=off).  */
  const char *no_tstvratio_help; /**< @brief If given fixed ts/tv ratios will not be used help description.  */
  float fixfactor_arg;	/**< @brief Float specifying what factor to use for saturated data. If not given -1 in the entry. (default='1').  */
  char * fixfactor_orig;	/**< @brief Float specifying what factor to use for saturated data. If not given -1 in the entry. original value given at command line.  */
  const char *fixfactor_help; /**< @brief Float specifying what factor to use for saturated data. If not given -1 in the entry. help description.  */
  int datasets_arg;	/**< @brief nr of datasets in input. Is used if the input format is clustalw (default='1').  */
  char * datasets_orig;	/**< @brief nr of datasets in input. Is used if the input format is clustalw original value given at command line.  */
  const char *datasets_help; /**< @brief nr of datasets in input. Is used if the input format is clustalw help description.  */
  int no_counts_flag;	/**< @brief Tree file should not contain the counts (default=off).  */
  const char *no_counts_help; /**< @brief Tree file should not contain the counts help description.  */
  enum enum_method *method_arg;	/**< @brief reconstruction method to apply.  */
  char ** method_orig;	/**< @brief reconstruction method to apply original value given at command line.  */
  int method_min; /**< @brief reconstruction method to apply's minimum occurreces */
  int method_max; /**< @brief reconstruction method to apply's maximum occurreces */
  const char *method_help; /**< @brief reconstruction method to apply help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int outfile_given ;	/**< @brief Whether outfile was given.  */
  unsigned int input_format_given ;	/**< @brief Whether input-format was given.  */
  unsigned int output_format_given ;	/**< @brief Whether output-format was given.  */
  unsigned int distance_function_given ;	/**< @brief Whether distance-function was given.  */
  unsigned int bootstraps_given ;	/**< @brief Whether bootstraps was given.  */
  unsigned int no_incl_orig_given ;	/**< @brief Whether no-incl-orig was given.  */
  unsigned int seed_given ;	/**< @brief Whether seed was given.  */
  unsigned int no_ambiguities_given ;	/**< @brief Whether no-ambiguities was given.  */
  unsigned int no_ambig_resolve_given ;	/**< @brief Whether no-ambig-resolve was given.  */
  unsigned int no_transprob_given ;	/**< @brief Whether no-transprob was given.  */
  unsigned int ambiguity_frequency_model_given ;	/**< @brief Whether ambiguity-frequency-model was given.  */
  unsigned int tstvratio_given ;	/**< @brief Whether tstvratio was given.  */
  unsigned int pyrtvratio_given ;	/**< @brief Whether pyrtvratio was given.  */
  unsigned int no_tstvratio_given ;	/**< @brief Whether no-tstvratio was given.  */
  unsigned int fixfactor_given ;	/**< @brief Whether fixfactor was given.  */
  unsigned int datasets_given ;	/**< @brief Whether datasets was given.  */
  unsigned int no_counts_given ;	/**< @brief Whether no-counts was given.  */
  unsigned int method_given ;	/**< @brief Whether method was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char * const *argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char * const *argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char * const *argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);

extern char *cmdline_parser_input_format_values[] ;	/**< @brief Possible values for input-format.  */
extern char *cmdline_parser_output_format_values[] ;	/**< @brief Possible values for output-format.  */
extern char *cmdline_parser_distance_function_values[] ;	/**< @brief Possible values for distance-function.  */
extern char *cmdline_parser_ambiguity_frequency_model_values[] ;	/**< @brief Possible values for ambiguity-frequency-model.  */
extern char *cmdline_parser_method_values[] ;	/**< @brief Possible values for method.  */


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* FASTDIST_GENGETOPT_H */
