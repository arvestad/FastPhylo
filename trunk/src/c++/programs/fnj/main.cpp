#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>

#include "config.h"
#include "file_utils.hpp"
#include "NeighborJoining.hpp"

#include "log_utils.hpp"
#include "fnj_gengetopt.h"

#include "DataInputStream.hpp"
#include "DataOutputStream.hpp"
#include "XmlOutputStream.hpp"
#include "Extrainfos.hpp"
#include "fileFormatSchema.hpp"

#ifdef WITH_LIBXML
#include "XmlInputStream.hpp"
#endif // WITH_LIBXML


using namespace std;


//
// Builds trees using the supplied distance methods.
// Each created tree is added to the tree2count map.
//
void
buildTrees(StrDblMatrix &dm, tree2int_map &tree2count, std::vector<NJ_method> &methods, str2int_hashmap &name2id){
  SequenceTree tree;
  for(size_t i=0 ; i<methods.size() ; i++){
    computeNJTree(dm,tree,methods[i]);
    tree.makeCanonical(name2id);
    tree2int_map::iterator iter = tree2count.find(tree);
    if(iter!=tree2count.end())
      (*iter).second++;
    else
      tree2count[tree] = 1;
  }
}

int
main(int argc,
     char **argv){

  gengetopt_args_info args_info;
  TRY_EXCEPTION();

  if (cmdline_parser (argc, argv, &args_info) != 0)
    exit(EXIT_FAILURE);

#ifndef WITH_LIBXML
  if ( args_info.input_format_arg == input_format_arg_xml ) {
     cerr << "The software was built with WITH_LIBXML=OFF. Please rebuild it if you want XML functionality." << endl; exit(EXIT_FAILURE);
  }
#endif // WITH_LIBXML

  if ( args_info.print_relaxng_input_given && args_info.print_relaxng_output_given ) {
     cerr << "error: --print-relaxng-input and --print-relaxng-output can not be used at the same time" << endl; exit(EXIT_FAILURE);
  }

  if ( args_info.print_relaxng_input_given ) {  cout << fastphylo_distance_matrix_xml_relaxngstr << std::endl;  exit(EXIT_SUCCESS);   };
  if ( args_info.print_relaxng_output_given ) {  cout << fastphylo_tree_count_xml_relaxngstr << std::endl;  exit(EXIT_SUCCESS);   };

  //----------------------------------------------
  // DISTANCE METHODS
  std::vector<NJ_method> methods;


  if( args_info.number_of_runs_given && args_info.input_format_arg == input_format_arg_xml )
    { cerr << "error: --number-of-runs can not be used together with input format xml." << endl; exit(EXIT_FAILURE);   }

  switch ( args_info.method_arg )
    { 
    case method_arg_NJ    :   methods.push_back(NJ); break;
    case method_arg_FNJ   :   methods.push_back(FNJ); break;
    case method_arg_BIONJ :   methods.push_back(BIONJ); break;
    default: cerr << "error: method chosen not available" << endl; exit(EXIT_FAILURE);
    }

  bool printCounts = args_info.print_counts_flag;

  try {

  char * inputfilename = 0;
  char * outputfilename = 0;

  DataInputStream *istream;
  DataOutputStream *ostream;

  switch( args_info.inputs_num )
    {  case 0: break; /* inputfilename will be null and indicate stdin as input */
    case 1: inputfilename =  args_info.inputs[0]; break;
    default: cerr << "Error: you can at most specify one input filename" << endl; exit(EXIT_FAILURE);
  }

  if( args_info.outfile_given )
    {  outputfilename = args_info.outfile_arg;  }

  switch ( args_info.input_format_arg )
    {
    case input_format_arg_phylip_dm: istream = new PhylipDmInputStream(inputfilename);  break;
#ifdef WITH_LIBXML
    case input_format_arg_xml: istream = new XmlInputStream(inputfilename); break;
#endif // WITH_LIBXML
    default: exit(EXIT_FAILURE);
}

  switch ( args_info.output_format_arg )
    {
    case output_format_arg_newick: ostream = new TreeTextOutputStream(outputfilename);  break;
    case output_format_arg_xml: ostream = new XmlOutputStream(outputfilename); break;
    default: exit(EXIT_FAILURE);
}
  
       // THE DATA WE WILL PROCESS
      std::vector<Sequence> seqs;
      std::vector<std::string> names;
      std::vector<DNA_b128_String> b128seqs;
      Extrainfos extrainfos;

      bool latestReadSuccessful = true;

      std::vector<std::string> speciesnames;

      readstatus status;
      int run = 0;
      status = END_OF_RUN;
      while (   status == END_OF_RUN   && ( args_info.input_format_arg == input_format_arg_xml || run < args_info.number_of_runs_arg )) {
        std::string runId("");
        run++;
        tree2int_map tree2count((size_t)( args_info.bootstraps_arg * 1.3));
        StrDblMatrix dm;
        str2int_hashmap name2id;
	int i;
        for ( i = 0 ; ( status == END_OF_RUN || status == DM_READ ) && ( i < args_info.dm_per_run_arg   || args_info.input_format_arg == input_format_arg_xml ) ; i++ ){

	  if (( status = istream->readDM( dm, names, runId, extrainfos)) != DM_READ ) { 
	    break;
	  };
     
    	  for(size_t namei=0 ; namei<dm.getSize() ; namei++ ) {
	     name2id[dm.getIdentifier(namei)] = namei;
	  }

	  buildTrees(dm, tree2count, methods,name2id);
        }
        if ( status == END_OF_RUN || i == args_info.dm_per_run_arg ) {
	  ostream->print(tree2count,printCounts, runId, names, extrainfos);
	}

      }//end run loop

      delete ostream;
      delete istream;
  }

  catch(...){
    throw;
  }

  CATCH_EXCEPTION();
  cmdline_parser_free(&args_info);
  return 0;
}

