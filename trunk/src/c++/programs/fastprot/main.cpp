#include "fastprot_gengetopt.h"
#include "log_utils.hpp"
#include "file_utils.hpp"
#include "PhylipMaInputStream.hpp"
#include "DataInputStream.hpp"
#include "FastaInputStream.hpp"
#include "BinaryDmOutputStream.hpp"
#include "DataOutputStream.hpp"
#include "PhylipDmOutputStream.hpp"
#include "XmlOutputStream.hpp"
#include "fileFormatSchema.hpp"
#include "ProtDistCalc.hpp"
#include "ProtSeqUtils.hpp"
#include "../../DistanceMatrix.hpp"
#include <string>
#include <vector>
#include <unistd.h>

#ifdef WITH_LIBXML
#include "XmlInputStream.hpp"
#endif //WITH_LIBXML

using namespace std;

int main (int argc, char **argv) {
  if(isatty(STDIN_FILENO) && argc==1) {
    cout<<"No input data or parameters. Use -h,--help for more information"<<endl;
    exit(EXIT_FAILURE);
  }
  gengetopt_args_info args_info;
  TRY_EXCEPTION();
  prot_sequence_translation_model trans_model;
#ifndef WITH_LIBXML
  if (args_info.input_format_arg == input_format_arg_xml){
    cerr << "The software was built with WITH_LIBXML=OFF. Please rebuild it if you want XML functionality." << endl;
    exit(EXIT_FAILURE);
  }
#endif // WITH_LIBXML
  if (cmdline_parser(argc, argv, &args_info) != 0)
    exit(EXIT_FAILURE);
  if (args_info.print_relaxng_input_given && args_info.print_relaxng_output_given) {
    cerr << "error: --print-relaxng-input and --print-relaxng-output can not be used at the same time" << endl;
    exit(EXIT_FAILURE);
  }
  if (args_info.print_relaxng_input_given) { 
    cout << fastphylo_prot_sequence_xml_relaxngstr << endl;
    exit(EXIT_SUCCESS);
  }
  if (args_info.print_relaxng_output_given) {
    cout << fastphylo_distance_matrix_xml_relaxngstr << endl;
    exit(EXIT_SUCCESS);
  }

//--------------------------------------------------------------
// Read translation model
// prot_sequence_translation_model trans_model;
  if (! args_info.model_file_given) {
    switch (args_info.distance_function_arg) {
      case distance_function_arg_ID:
        trans_model.model = id;
        break;
      case distance_function_arg_JC:
        trans_model.model = jc;
        break;
      case distance_function_arg_JCK:
        trans_model.model = jck;
        break;
      case distance_function_arg_JCSS:
        trans_model.model = jcss;
        break;
      case distance_function_arg_WAG:
        trans_model.model = wag;
        break;
      case distance_function_arg_JTT:
        trans_model.model = jtt;
        break;
      case distance_function_arg_DAY:
        trans_model.model = day;
        break;
      case distance_function_arg_ARVE:
        trans_model.model = arve;
        break;
      case distance_function_arg_MVR:
        trans_model.model = mvr;
        break;
      case distance_function_arg_LG:
              trans_model.model = lg;
              break;
      default:
        cerr << "error: model chosen not available" << endl;
        exit(EXIT_FAILURE);
    }
  } else {
    // read file from args_info.model_file_arg
  }
  if (args_info.maximum_likelihood_given && 
      (trans_model.model == id || trans_model.model == jc ||
      trans_model.model == jck || trans_model.model == jcss ||
      args_info.sd_given)) {
    cerr << "error: --maximum-likelihood can not be used with --distance-function=ID, JC, JCK or JCSS or --sd" << endl;
    exit(EXIT_FAILURE);
  }
  if (args_info.sd_given && 
      (trans_model.model == id || trans_model.model == jc ||
      trans_model.model == jck || trans_model.model == jcss ||
      args_info.maximum_likelihood_given)) {
    cerr << "error: --sd can not be used with --distance-function=ID, JC, JCK or JCSS or --maximum-likelihood" << endl;
    exit(EXIT_FAILURE);
  }
  if (args_info.sd_given) {
    // Remove this two lines when sd is working again
    cerr << "Aborting: The std dev feature is not yet implemented." << endl;
    exit(1);
  }
  trans_model.step_size = args_info.speed_arg;
  if (args_info.pfam_given)
    trans_model.tp = norm;
  else 
    trans_model.tp = flat;
  trans_model.sd = args_info.sd_given;
  trans_model.ml = args_info.maximum_likelihood_given;
  bool remove_indels = args_info.remove_indels_given;
  int ndatasets= 1;
//----------------------------------------------
// BOOTSTRAPPING
  int numboot = args_info.bootstraps_arg;
  bool no_incl_orig = args_info.no_incl_orig_given;
  if (numboot==0 && no_incl_orig) {
    cerr << "error: No output. No bootstrap or original data will be written" << endl;
    exit(EXIT_FAILURE);
    }
  if ( args_info.seed_given ) 
    srand((unsigned int )args_info.seed_arg);
  else
    srand((unsigned int)time(NULL));
  try {
    char *inputfilename = NULL;
    char *outputfilename = NULL;
    DataInputStream *istream;
    DataOutputStream *ostream;
    switch( args_info.inputs_num ) {
      case 0:
        break; /* inputfilename will be null and indicate stdin as input */
      case 1:
        inputfilename =  args_info.inputs[0];
        break;
      default:
        cerr << "Error: you can at most specify one input filename" << endl;
        exit(EXIT_FAILURE);
    }
    if( args_info.outfile_given )
      outputfilename = args_info.outfile_arg;
    switch (args_info.input_format_arg) {
      case input_format_arg_fasta:
        istream = new FastaInputStream(inputfilename);
        break;
      case input_format_arg_phylip:
        istream = new PhylipMaInputStream(inputfilename);
        break;
#ifdef WITH_LIBXML
      case input_format_arg_xml:
        istream = new XmlInputStream(inputfilename);
        break;
#endif // WITH_LIBXML
      default:
        exit(EXIT_FAILURE);
    }
    bool binary_format_type=args_info.memory_efficient_given;
    switch (args_info.output_format_arg) {
      case output_format_arg_phylip:
        ostream = new PhylipDmOutputStream(outputfilename);
        break;
      case output_format_arg_xml:
        ostream = new XmlOutputStream(outputfilename);
        break;
      case output_format_arg_binary:
        ostream = new BinaryDmOutputStream(outputfilename);
        binary_format_type=true;
        break;
      default:
        exit(EXIT_FAILURE);
    }

		StrDblMatrix dm;
		StrDblMatrix sdm;
		vector<Sequence> seqs;
		vector<string> names;
		Extrainfos extrainfos;

		for ( int ds = 0 ; ds < ndatasets || args_info.input_format_arg == input_format_arg_xml ; ds++ ){
			string runId("");
			if (! istream->read(seqs, runId, names, extrainfos))
				break;
			if (remove_indels)
				remove_gaps(seqs);
			if (!no_incl_orig) {
				if (trans_model.sd)
					calculate_distances(seqs, dm, trans_model, sdm);
				else
					calculate_distances(seqs, dm, trans_model);
				dm.setIdentifiers(names);
				ostream->printStartRun(names, runId, extrainfos);
				ostream->printHeader(seqs.size());
				ostream->print(dm);
				if (trans_model.sd && !binary_format_type)
					ostream->printSD(sdm);
			}
			// Bootstrapping
			for (int b=0; b < numboot; b++) {
				vector<Sequence> bseqs;
				bootstrap_sequences(seqs, bseqs);
				if (trans_model.sd)
					calculate_distances(bseqs, dm, trans_model, sdm);
				else
					calculate_distances(bseqs, dm, trans_model);
				dm.setIdentifiers(names);
				ostream->print(dm);
				if (trans_model.sd && !binary_format_type)
					ostream->printSD(sdm);
			}
			if (!binary_format_type)
				ostream->printEndRun();
		}
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
