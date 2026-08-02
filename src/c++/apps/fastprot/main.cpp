#include "config.h"
#include "fastprot_gengetopt.h"
#include "fastphylo/core/log_utils.hpp"
#include "fastphylo/core/file_utils.hpp"
#include "PhylipMaInputStream.hpp"
#include "DataInputStream.hpp"
#include "FastaInputStream.hpp"
#include "fastphylo/io/BinaryDmOutputStream.hpp"
#include "fastphylo/io/DataOutputStream.hpp"
#include "fastphylo/io/PhylipDmOutputStream.hpp"
#include "fastphylo/io/XmlOutputStream.hpp"
#include "fastphylo/io/fileFormatSchema.hpp"
#include "fastphylo/protein/ProtDistCalc.hpp"
#include "fastphylo/protein/ProtSeqUtils.hpp"
#include "fastphylo/core/DistanceMatrix.hpp"
#include <memory>
#include <string>
#include <vector>
#include <unistd.h>

#ifdef WITH_LIBXML
#include "XmlInputStream.hpp"
#endif //WITH_LIBXML

using namespace std;

// Lint Phase 5 (lint_phase5_refactors.md): main() was a single 220-line
// function (cognitive complexity 59, threshold 25). Extracted below, in
// the order main() calls them.
static void handleEarlyExitFlags(const gengetopt_args_info &args_info) {
#ifndef WITH_LIBXML
  if (args_info.input_format_arg == input_format_arg_xml){
    cerr << "The software was built with WITH_LIBXML=OFF. Please rebuild it if you want XML functionality." << endl;
    exit(EXIT_FAILURE);
  }
#endif // WITH_LIBXML
  if ((args_info.print_relaxng_input_given != 0u) && (args_info.print_relaxng_output_given != 0u)) {
    cerr << "error: --print-relaxng-input and --print-relaxng-output can not be used at the same time" << endl;
    exit(EXIT_FAILURE);
  }
  if (args_info.print_relaxng_input_given != 0u) {
    cout << fastphylo_prot_sequence_xml_relaxngstr << endl;
    exit(EXIT_SUCCESS);
  }
  if (args_info.print_relaxng_output_given != 0u) {
    cout << fastphylo_distance_matrix_xml_relaxngstr << endl;
    exit(EXIT_SUCCESS);
  }
}

static prot_sequence_translation_model buildTranslationModel(const gengetopt_args_info &args_info) {
  prot_sequence_translation_model trans_model;
  if (args_info.model_file_given == 0u) {
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
  if ((args_info.maximum_likelihood_given != 0u) &&
      (trans_model.model == id || trans_model.model == jc ||
      trans_model.model == jck || trans_model.model == jcss ||
      (args_info.sd_given != 0u))) {
    cerr << "error: --maximum-likelihood can not be used with --distance-function=ID, JC, JCK or JCSS or --sd" << endl;
    exit(EXIT_FAILURE);
  }
  if ((args_info.sd_given != 0u) &&
      (trans_model.model == id || trans_model.model == jc ||
      trans_model.model == jck || trans_model.model == jcss ||
      (args_info.maximum_likelihood_given != 0u))) {
    cerr << "error: --sd can not be used with --distance-function=ID, JC, JCK or JCSS or --maximum-likelihood" << endl;
    exit(EXIT_FAILURE);
  }
  if (args_info.sd_given != 0u) {
    // Remove this two lines when sd is working again
    cerr << "Aborting: The std dev feature is not yet implemented." << endl;
    exit(1);
  }
  trans_model.step_size = args_info.speed_arg;
  if (args_info.pfam_given != 0u) {
    trans_model.tp = norm;
  } else {
    trans_model.tp = flat;
  }
  trans_model.sd = (args_info.sd_given != 0u);
  trans_model.ml = (args_info.maximum_likelihood_given != 0u);
  return trans_model;
}

static std::unique_ptr<DataInputStream> buildInputStream(const gengetopt_args_info &args_info) {
  char *inputfilename = nullptr;
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
  switch (args_info.input_format_arg) {
    case input_format_arg_fasta:
      return std::make_unique<FastaInputStream>(inputfilename);
    case input_format_arg_phylip:
      return std::make_unique<PhylipMaInputStream>(inputfilename);
#ifdef WITH_LIBXML
    case input_format_arg_xml:
      return std::make_unique<XmlInputStream>(inputfilename);
#endif // WITH_LIBXML
    default:
      exit(EXIT_FAILURE);
  }
}

static std::unique_ptr<DataOutputStream> buildOutputStream(const gengetopt_args_info &args_info) {
  char *outputfilename = nullptr;
  if( args_info.outfile_given != 0u ) {
    outputfilename = args_info.outfile_arg;
  }
  switch (args_info.output_format_arg) {
    case output_format_arg_phylip:
      return std::make_unique<PhylipDmOutputStream>(outputfilename);
    case output_format_arg_xml:
      return std::make_unique<XmlOutputStream>(outputfilename);
    case output_format_arg_binary:
      return std::make_unique<BinaryDmOutputStream>(outputfilename);
    default:
      exit(EXIT_FAILURE);
  }
}

// Shared by processRuns()'s original-data and bootstrap-replicate
// cases below: pick sd/non-sd calculate_distances(), then stamp
// identifiers - everything else around it (whether printStartRun()/
// printHeader() run first) differs between the two call sites, so
// only this common part is factored out.
static void calculateDistances(prot_sequence_translation_model &trans_model, vector<Sequence> &seqs,
                                StrDblMatrix &dm, StrDblMatrix &sdm, vector<string> &names) {
  if (trans_model.sd) {
    calculate_distances(seqs, dm, trans_model, sdm);
  } else {
    calculate_distances(seqs, dm, trans_model);
  }
  dm.setIdentifiers(names);
}

static void printDistances(DataOutputStream &ostream, prot_sequence_translation_model &trans_model,
                            StrDblMatrix &dm, StrDblMatrix &sdm, bool binary_format_type) {
  ostream.print(dm);
  if (trans_model.sd && !binary_format_type) {
    ostream.printSD(sdm);
  }
}

// The dataset/bootstrap run loop - reads each dataset, computes and
// prints its distance matrix, then numboot bootstrap replicates of it.
static void processRuns(DataInputStream &istream, DataOutputStream &ostream,
                         prot_sequence_translation_model &trans_model, int ndatasets,
                         int numboot, bool no_incl_orig, bool remove_indels,
                         bool binary_format_type, int input_format_arg) {
  StrDblMatrix dm;
  StrDblMatrix sdm;
  vector<Sequence> seqs;
  vector<string> names;
  Extrainfos extrainfos;

  for ( int ds = 0 ; ds < ndatasets || input_format_arg == input_format_arg_xml ; ds++ ){
    string runId;
    if (! istream.read(seqs, runId, names, extrainfos)) {
      break;
    }
    if (remove_indels) {
      remove_gaps(seqs);
    }
    if (!no_incl_orig) {
      calculateDistances(trans_model, seqs, dm, sdm, names);
      ostream.printStartRun(names, runId, extrainfos);
      ostream.printHeader(seqs.size());
      printDistances(ostream, trans_model, dm, sdm, binary_format_type);
    }
    // Bootstrapping
    for (int b=0; b < numboot; b++) {
      vector<Sequence> bseqs;
      bootstrap_sequences(seqs, bseqs);
      calculateDistances(trans_model, bseqs, dm, sdm, names);
      printDistances(ostream, trans_model, dm, sdm, binary_format_type);
    }
    if (!binary_format_type) {
      ostream.printEndRun();
    }
  }
}

int main (int argc, char **argv) {
  if((isatty(STDIN_FILENO) != 0) && argc==1) {
    cout<<"No input data or parameters. Use -h,--help for more information"<<endl;
    exit(EXIT_FAILURE);
  }
  gengetopt_args_info args_info;
  TRY_EXCEPTION();
  if (cmdline_parser(argc, argv, &args_info) != 0) {
    exit(EXIT_FAILURE);
}
  handleEarlyExitFlags(args_info);

  prot_sequence_translation_model trans_model = buildTranslationModel(args_info);
  bool remove_indels = args_info.remove_indels_given != 0u;
  int ndatasets= 1;
//----------------------------------------------
// BOOTSTRAPPING
  int numboot = args_info.bootstraps_arg;
  bool no_incl_orig = args_info.no_incl_orig_given != 0u;
  if (numboot==0 && no_incl_orig) {
    cerr << "error: No output. No bootstrap or original data will be written" << endl;
    exit(EXIT_FAILURE);
    }
  if ( args_info.seed_given != 0u ) {
    srand(static_cast<unsigned int>(args_info.seed_arg));
  } else {
    // NOLINTNEXTLINE(bugprone-random-generator-seed) - bootstrap resampling, not a security context; time-seeding when no explicit --seed is the intended default.
    srand(static_cast<unsigned int>(time(nullptr)));
}
  try {
    std::unique_ptr<DataInputStream> istream = buildInputStream(args_info);
    std::unique_ptr<DataOutputStream> ostream = buildOutputStream(args_info);
    bool binary_format_type = (args_info.memory_efficient_given != 0u) ||
                               (args_info.output_format_arg == output_format_arg_binary);

    processRuns(*istream, *ostream, trans_model, ndatasets, numboot, no_incl_orig,
                remove_indels, binary_format_type, args_info.input_format_arg);
  }
  catch(...){
    throw;
  }
  CATCH_EXCEPTION()
  catch(...){
    std::cerr << "Unknown (non-Exception) error" << std::endl;
  }
  cmdline_parser_free(&args_info);
  return 0;
}
