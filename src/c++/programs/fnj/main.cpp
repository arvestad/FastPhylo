

#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>

#include "file_utils.hpp"
#include "NeighborJoining.hpp"



#include "log_utils.hpp"
#include "fnj_gengetopt.h"
#include "NeighborJoining.hpp"


#include "DataInputStream.hpp"
#include "XmlInputStream.hpp"
#include "DataOutputStream.hpp"
#include "XmlOutputStream.hpp"

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

  TRY_EXCEPTION();

  gengetopt_args_info args_info;

  if (cmdline_parser (argc, argv, &args_info) != 0)
    exit(EXIT_FAILURE);

  if ( args_info.print_relaxng_given ) {  cout << relaxngstr << std::endl;  exit(EXIT_SUCCESS);   };

  //----------------------------------------------
  // DISTANCE METHODS
  std::vector<NJ_method> methods;


  if( args_info.number_of_runs_given && args_info.input_format_arg == input_format_arg_xml )
    { cerr << "error: --number-of-runs can not be used together with input format xml." << endl; exit(EXIT_FAILURE);   }

  if( ! args_info.method_given ) 
    { cerr << "error: at least one --method must be given" << endl; exit(EXIT_FAILURE);   }




  for (int i = 0; i < args_info.method_given; ++i)
    {
      switch ( args_info.method_arg[i] )
	{ 
	case method_arg_NJ    :   methods.push_back(NJ); break;
	case method_arg_FNJ   :   methods.push_back(FNJ); break;
	case method_arg_BIONJ :   methods.push_back(BIONJ); break;
	default: cerr << "error: method chosen not available" << endl; exit(EXIT_FAILURE);
	}
    }

  //--------------
  bool noCounts = args_info.no_counts_flag;



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
    case input_format_arg_phylip_dm: istream = new PhylipMaInputStream(inputfilename);  break;
    case input_format_arg_xml: istream = new XmlInputStream(inputfilename); break;
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
      std::vector<string> names;
      std::vector<DNA_b128_String> b128seqs;

      //for each dataset in the files

      //      for ( int ds = 0 ; ds < ndatasets || args_info.input_format_arg == input_format_arg_xml ; ds++ ){

      bool latestReadSuccessful = true;

      std::vector<std::string> speciesnames;

      readstatus status;
      int run = 0;
      status = END_OF_RUN;
      while (   status == END_OF_RUN   && ( args_info.input_format_arg == input_format_arg_xml || run < args_info.number_of_runs_arg )) {
        run++;
        tree2int_map tree2count((size_t)( args_info.bootstraps_arg * 1.3));
        StrDblMatrix dm;
        str2int_hashmap name2id;

        for ( int i = 0 ; ( status == END_OF_RUN || status == DM_READ ) && ( i < args_info.dm_per_run_arg - 1 || args_info.input_format_arg == input_format_arg_xml ) ; i++ ){

	  //	   applyFixFactor(dm,fixfactor);
	  if ( ( status = istream->readDM( dm )) != DM_READ ) { 
	    break;
	  };
     
	for(size_t namei=0 ; namei<dm.getSize() ; namei++ )
	  name2id[dm.getIdentifier(namei)] = namei;

	  buildTrees(dm, tree2count, methods,name2id);
        }

        ostream->print(tree2count,noCounts);

      }//end run loop

      delete ostream;
      delete istream;
  }

  catch(...){
    throw;
  }

  CATCH_EXCEPTION();
  return 0;
}

