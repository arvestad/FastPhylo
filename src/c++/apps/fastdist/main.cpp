///////////////////////////////////////////////
//                                           //
// File: fastdist.cpp                        //
//                                           //
// Created on: <14-Mar-2005 14:21:14 cherub> //
//                                           //
// Author: Isaac Elias,Mehmood Alam Khan     //
// Email: isaac@nada.kth.se, malagori@kth.se //
///////////////////////////////////////////////

#include "fastphylo/dna/Sequences2DistanceMatrix.hpp"

#include <memory>
#include <string>
#include <iostream>
#include <ctime>
#include <fstream>
#include <cassert>
#include <map>
// Windows' CRT defaults stdout to text mode, translating every '\n'
// this project's own code writes into "\r\n" - this project's output
// is meant to be byte-identical across platforms (RunExamples.sh
// byte-diffs it), and the code always writes '\n' explicitly, never
// relies on CRT translation, so binary mode is the correct default
// here, not just a workaround. _setmode()/_fileno() are declared in
// <io.h>/<fcntl.h>, MSVC-only.
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
static void setStdoutBinaryMode() { _setmode(_fileno(stdout), _O_BINARY); }
#else
static void setStdoutBinaryMode() {}
#endif

#include "config.h"
#include "fastphylo/core/file_utils.hpp"
#include <iomanip>
#include "fastphylo/core/log_utils.hpp"
#include "fastphylo/io/BinaryDmOutputStream.hpp"
#include "CliHelpFormatter.hpp"
#include "fastphylo/dna/NeighborJoining.hpp"
#include "DataInputStream.hpp"
#include "PhylipMaInputStream.hpp"
#include "FastaInputStream.hpp"
#include "fastphylo/io/DataOutputStream.hpp"
#include "fastphylo/io/Extrainfos.hpp"
#include "fastphylo/io/fileFormatSchema.hpp"
#include "fastphylo/io/XmlOutputStream.hpp"
#include "fastphylo/io/PhylipDmOutputStream.hpp"

#ifdef WITH_LIBXML
#include "XmlInputStream.hpp"
#endif // WITH_LIBXML

using namespace std;

// gengetopt_migration_plan.md, Phase C: fastdist is the second app
// migrated off gengetopt onto CLI11 (fnj was first, Phase B).
// FastdistOptions replaces gengetopt_args_info, same shape as fnj's
// FnjOptions - values plus "_given" bools for the options whose
// *presence*, not just their value, changes behavior (fixfactor is
// the interesting one: it has a default of 1.0 either way, but
// useFixFactor below only applies it when --fixfactor was actually
// given). --distance-function binds directly to
// dna_pairwise_sequence_likelihood.hpp's existing sequence_model
// domain enum via CLI::CheckedTransformer, same trick as fnj's
// --method/NJ_method. input-format/output-format/ambiguity-frequency
// -model get their own small scoped enums (fastdist's value sets
// don't match fnj's input-format/output-format enums, so those aren't
// reused).
enum class InputFormat { Fasta, Phylip, Xml };
enum class OutputFormat { Phylip, Xml, Binary };
enum class AmbiguityFreqModel { Uni, Base };

// The .ggo's own "version" line (independently pinned there, not tied
// to the top-level CMakeLists.txt's PACKAGE_VERSION="1.0.3") - kept
// exactly as-is, same reasoning as fnj's FNJ_VERSION.
static constexpr const char *FASTDIST_VERSION = "fastdist 1.0.10";

struct FastdistOptions {
	std::string outfile;
	bool outfile_given = false;
	InputFormat input_format = InputFormat::Fasta;
	bool memory_efficient = false;
	OutputFormat output_format = OutputFormat::Phylip;
	sequence_model distance_function = K2P;
	int bootstraps = 0;
	bool no_incl_orig = false;
	int seed = 0;
	bool seed_given = false;
	bool no_ambiguities = false;
	bool no_ambig_resolve = false;
	bool no_transprob = false;
	AmbiguityFreqModel ambiguity_frequency_model = AmbiguityFreqModel::Uni;
	float tstvratio = 2.0F;
	float pyrtvratio = 2.0F;
	bool no_tstvratio = false;
	float fixfactor = 1.0F;
	bool fixfactor_given = false;
	int number_of_runs = 1;
	bool number_of_runs_given = false;
	bool print_relaxng_input = false;
	bool print_relaxng_output = false;
	std::string inputfile;
	bool inputfile_given = false;
};

// Builds the CLI11 App and parses argv into a FastdistOptions. Exits
// the process directly on --help/--version/a parse error, same
// approach as fnj's parseArgs() (Phase B).
static FastdistOptions parseArgs(int argc, char **argv) {
	FastdistOptions opts;
	CLI::App app{"Computes distance matrices out of multialignments", "fastdist"};
	app.formatter(std::make_shared<FastphyloHelpFormatter>());
	app.footer("\nIf FILE is not specified the input is read from stdin\n\n"
	           "Example usage of this program can be found at its home page\n"
	           "http://fastphylo.sourceforge.net/\n");
	app.set_version_flag("-V,--version", FASTDIST_VERSION);

	app.add_option("-o,--outfile", opts.outfile,
	                "output filename. If not specifed, output is written to stdout")
	    ->type_name("filename");

	static const std::map<std::string, InputFormat> input_format_map{
	    {"fasta", InputFormat::Fasta}, {"phylip", InputFormat::Phylip}, {"xml", InputFormat::Xml}};
	app.add_option("-I,--input-format", opts.input_format,
	                "input format. xml means the Fastphylo sequence XML format "
	                "(possible values: fasta, phylip, xml; default: fasta)")
	    ->transform(CLI::CheckedTransformer(input_format_map).description(""));

	app.add_flag("-e,--memory-efficient", opts.memory_efficient,
	             "memory efficient. Use less memory space and fast implementation. "
	             "Only used with fasta and phylip format");

	static const std::map<std::string, OutputFormat> output_format_map{
	    {"phylip", OutputFormat::Phylip}, {"xml", OutputFormat::Xml}, {"binary", OutputFormat::Binary}};
	app.add_option("-O,--output-format", opts.output_format,
	                "output format. xml means the Fastphylo distance matrix XML format "
	                "(possible values: phylip, xml, binary; default: phylip)")
	    ->transform(CLI::CheckedTransformer(output_format_map).description(""));

	static const std::map<std::string, sequence_model> distance_function_map{
	    {"JC", JC}, {"K2P", K2P}, {"TN93", TN93}, {"HAMMING", HAMMING_DISTANCE}};
	app.add_option("-D,--distance-function", opts.distance_function,
	                "Distance function (possible values: JC, K2P, TN93, HAMMING; default: K2P)")
	    ->transform(CLI::CheckedTransformer(distance_function_map).description(""));

	app.add_option("-b,--bootstraps", opts.bootstraps, "Bootstrap num times and create matrix for each");

	app.add_flag("-k,--no-incl-orig", opts.no_incl_orig,
	             "If the distance matrix from the original sequences should not be included");

	CLI::Option *seed_opt = app.add_option(
	    "-s,--seed", opts.seed, "Random seed. If not specified the current timestamp will be used");

	app.add_flag("-A,--no-ambiguities", opts.no_ambiguities, "Ignore ambiguities");

	app.add_flag("-R,--no-ambig-resolve", opts.no_ambig_resolve,
	             "Specifies that ambigious symbols should not be resolved by nearest neighbor");

	app.add_flag("-t,--no-transprob", opts.no_transprob,
	             "Specifies that the transition probabilities should not be used in the ambiguity model");

	static const std::map<std::string, AmbiguityFreqModel> ambiguity_freq_map{
	    {"UNI", AmbiguityFreqModel::Uni}, {"BASE", AmbiguityFreqModel::Base}};
	app.add_option("-a,--ambiguity-frequency-model", opts.ambiguity_frequency_model,
	                "Ambiguity frequency model (possible values: UNI, BASE; default: UNI)")
	    ->transform(CLI::CheckedTransformer(ambiguity_freq_map).description(""));

	app.add_option("-T,--tstvratio", opts.tstvratio,
	               "Transition/transvertion ratio for purine transitions (for the TN model)");

	app.add_option("-P,--pyrtvratio", opts.pyrtvratio,
	               "Transition/transvertion ratio for  pyrimidines transitions (for the TN model)");

	app.add_flag("-N,--no-tstvratio", opts.no_tstvratio, "If given fixed ts/tv ratios will not be used");

	CLI::Option *fixfactor_opt = app.add_option(
	    "-F,--fixfactor", opts.fixfactor,
	    "Float specifying what factor to use for saturated data. If not given -1 in the entry.");

	CLI::Option *runs_opt = app.add_option(
	    "-r,--number-of-runs", opts.number_of_runs,
	    "nr of runs (datasets) in input. This option is only used if the input format is "
	    "phylip_multialignment.");

	app.add_flag("-p,--print-relaxng-input", opts.print_relaxng_input,
	             "print the Relax NG schema for the XML input format (Fastphylo sequence XML "
	             "format) and then exit");

	app.add_flag("-w,--print-relaxng-output", opts.print_relaxng_output,
	             "print the Relax NG schema for the XML output format (Fastphylo distance "
	             "matrix XML format) and then exit.");

	CLI::Option *input_opt =
	    app.add_option("FILE", opts.inputfile, "input file. If not specified, input is read from stdin");

	try {
		app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
		app.exit(e);
		// gengetopt exits 1 uniformly on any parse problem; CLI11's own
		// exit codes vary per error type. --help/--version (CLI::Success
		// and its subclasses) are the one case that's actually 0.
		std::exit(dynamic_cast<const CLI::Success *>(&e) != nullptr ? EXIT_SUCCESS : EXIT_FAILURE);
	}

	opts.outfile_given = app["--outfile"]->count() > 0;
	opts.seed_given = seed_opt->count() > 0;
	opts.fixfactor_given = fixfactor_opt->count() > 0;
	opts.number_of_runs_given = runs_opt->count() > 0;
	opts.inputfile_given = input_opt->count() > 0;
	return opts;
}

// Lint Phase 5 (lint_phase5_refactors.md): main() was a single 394-line
// function (cognitive complexity 181, threshold 25). Extracted below,
// in the order main() calls them.

// Handles --print-relaxng-input/--print-relaxng-output/--number-of-runs'
// argument-combination checks and the two exit(EXIT_SUCCESS) relaxng-dump
// paths - none of this depends on anything computed later in main().
static void handleEarlyExitFlags(const FastdistOptions &opts) {
#ifndef WITH_LIBXML
	if ( opts.input_format == InputFormat::Xml ) {
		cerr << "The software was built with WITH_LIBXML=OFF. Please rebuild it if you want XML functionality." << endl; exit(EXIT_FAILURE);
	}
#endif // WITH_LIBXML

	if ( opts.print_relaxng_input && opts.print_relaxng_output ) {
		cerr << "error: --print-relaxng-input and --print-relaxng-output can not be used at the same time" << endl; exit(EXIT_FAILURE);
	}

	if ( opts.print_relaxng_input ) {  cout << fastphylo_sequence_xml_relaxngstr << std::endl;  exit(EXIT_SUCCESS);   };
	if ( opts.print_relaxng_output ) {  cout << fastphylo_distance_matrix_xml_relaxngstr << std::endl;  exit(EXIT_SUCCESS);   };

	if ( opts.number_of_runs_given && opts.input_format != InputFormat::Phylip ) {
		cerr << "error: --number-of-runs can only be used together with --input-format=phylip " << endl; exit(EXIT_FAILURE);
	}
}

// Builds the sequence_translation_model from opts: evolutionary
// model choice, ts/tv ratios, and the ambiguity-handling flags.
static sequence_translation_model buildTranslationModel(const FastdistOptions &opts) {
	sequence_translation_model trans_model;

	trans_model.model = opts.distance_function;

	trans_model.no_tstvratio = opts.no_tstvratio;
	trans_model.tstvratio = opts.tstvratio;
	trans_model.pyrtvratio = opts.pyrtvratio;

	trans_model.no_ambiguities = opts.no_ambiguities;
	trans_model.no_ambig_resolve = opts.no_ambig_resolve;
	trans_model.no_transition_probs = opts.no_transprob;

	trans_model.use_base_freqs = (opts.ambiguity_frequency_model == AmbiguityFreqModel::Base);

	return trans_model;
}

// srand()'s seed: an explicit --seed for reproducibility, or the current
// time by default (bootstrap resampling, not a security context).
static void seedRandom(const FastdistOptions &opts) {
	if ( opts.seed_given ) {
		srand(static_cast<unsigned int>(opts.seed));
	} else {
		// NOLINTNEXTLINE(bugprone-random-generator-seed) - bootstrap resampling, not a security context; time-seeding when no explicit --seed is the intended default.
		srand(static_cast<unsigned int>(time(nullptr)));
	}
}

static std::unique_ptr<DataInputStream> buildInputStream(const FastdistOptions &opts) {
	char *inputfilename = nullptr;
	if ( opts.inputfile_given ) {
		inputfilename = const_cast<char *>(opts.inputfile.c_str());
	}

	switch ( opts.input_format )
	{
	case InputFormat::Fasta: return std::make_unique<FastaInputStream>(inputfilename);
	case InputFormat::Phylip : return std::make_unique<PhylipMaInputStream>(inputfilename);
#ifdef WITH_LIBXML
	case InputFormat::Xml: return std::make_unique<XmlInputStream>(inputfilename);
#endif // WITH_LIBXML
	default: exit(EXIT_FAILURE);
	}
}

static std::unique_ptr<DataOutputStream> buildOutputStream(const FastdistOptions &opts) {
	char *outputfilename = nullptr;
	if( opts.outfile_given )
	{  outputfilename = const_cast<char *>(opts.outfile.c_str());  }

	switch ( opts.output_format )
	{
	case OutputFormat::Phylip: return std::make_unique<PhylipDmOutputStream>(outputfilename);
	case OutputFormat::Xml: return std::make_unique<XmlOutputStream>(outputfilename);
	//Mehmood's Changes here : email: malagori@kth.se
	// binary_dm_format_plan.md: matricesPerRun is a single value fixed
	// for the whole invocation - numboot/no_incl_orig never change per
	// dataset even across ndatasets=opts.number_of_runs iterations.
	case OutputFormat::Binary: return std::make_unique<BinaryDmOutputStream>(
	    outputfilename, (opts.no_incl_orig ? 0 : 1) + static_cast<size_t>(opts.bootstraps));
	default: exit(EXIT_FAILURE);
	}
}

// Row-streaming output has two modes, distinguished at every call site
// by name instead of a bare true/false: Regular streams rows for
// --output-format=binary, MemoryEfficient additionally reuses buffers
// for --memory-efficient. Converted to the bool fillMatrixRow()/
// printRow() still expect only at the point of calling them.
enum class RowOutputMode { Regular, MemoryEfficient };

// Fills and prints every row of dm - the innermost loop of
// runRowStreamingOutput() below, which the original had hand-duplicated
// three times (no-bootstrap case, bootstrap's original-sequences case,
// bootstrap's per-replicate case) with only dm/dm2 differing.
static void fillAndPrintAllRows(DataOutputStream &ostream, StrFloRow &dm,
                                 std::vector<DNA_b128_String> &b128seqs,
                                 sequence_translation_model &trans_model,
                                 std::vector<string> &names, size_t numberOfSequences,
                                 bool useFixFactor, float fixfactor, RowOutputMode mode) {
	const bool mem_eff_flag = (mode == RowOutputMode::MemoryEfficient);
	for(size_t i = 0; i < numberOfSequences; ++i){
		fillMatrixRow(dm, b128seqs, trans_model, i, mem_eff_flag);
		dm.setIdentifier(names.at(i));
		if(useFixFactor) {
			applyFixFactorRow(dm,fixfactor);
		}
		ostream.printRow(dm, names.at(i), i, mem_eff_flag);
	}
}

// Row-streaming output path (StrFloRow, fillMatrixRow()/printRow()) -
// shared by the --output-format=binary and --memory-efficient branches,
// which turned out to be identical apart from mem_eff_flag once compared
// side by side (both originally hand-duplicated the same ~75 lines).
static void runRowStreamingOutput(DataInputStream &istream, DataOutputStream &ostream,
                                   sequence_translation_model &trans_model,
                                   int ndatasets, int numboot, bool no_incl_orig,
                                   bool useFixFactor, float fixfactor,
                                   InputFormat input_format, RowOutputMode mode) {
	StrFloRow dm;
	//open infile

	// THE DATA WE WILL PROCESS
	std::vector<Sequence> seqs;
	std::vector<string> names;
	std::vector<DNA_b128_String> b128seqs;

	Extrainfos extrainfos;

	//for each dataset in the files

	for ( int ds = 0 ; ds < ndatasets || input_format == InputFormat::Xml ; ds++ ){
		//no bootstrapping
		std::string runId;
		if ( !no_incl_orig && numboot == 0){//only need to create one distance matrix
			if ( ! istream.read(b128seqs,runId,names,extrainfos)) {
				break;
			}

			const size_t numberOfSequences = b128seqs.size();
			ostream.printStartRun(names,runId,extrainfos);
			ostream.printHeader(numberOfSequences);
			fillAndPrintAllRows(ostream, dm, b128seqs, trans_model, names, numberOfSequences, useFixFactor, fixfactor, mode);
		}
		//bootstrapping
		else{
			//TODO: implement bootstraping with the improved version of FastDist
			StrFloRow dm2;
			//read original sequences
			if ( ! istream.readSequences(seqs,runId,extrainfos)) { break;
}
			names.clear();names.reserve(seqs.size());
			for(auto & seq : seqs)
			{
				names.push_back(seq.name);
			}
			const size_t numberOfSequences = seqs.size();
			ostream.printStartRun(names,runId,extrainfos);
			ostream.printHeader(numberOfSequences);
			if ( !no_incl_orig ){//create the distance matrix for the original sequences
				Sequences2DNA_b128(seqs,b128seqs);
				fillAndPrintAllRows(ostream, dm2, b128seqs, trans_model, names, numberOfSequences, useFixFactor, fixfactor, mode);
			}
			//start the bootstrapping
			//	  vector<Sequence> bootsequences;
			for ( int b = 0 ; b < numboot ; b++ ){
				bootstrapSequences(seqs,b128seqs);
				std::cout << numberOfSequences << std::endl;
				ostream.printBootstrapSpliter(numberOfSequences);
				fillAndPrintAllRows(ostream, dm2, b128seqs, trans_model, names, numberOfSequences, useFixFactor, fixfactor, mode);
			}
		}
		ostream.printEndRun();
	}//end data set loop
}

// Sets dm's identifiers, applies the fix factor if requested, and prints
// it - runFullMatrixOutput()'s equivalent of fillAndPrintAllRows() above,
// hand-duplicated three times in the original for the same three cases.
// fillMatrix() itself stays a separate call at each site rather than
// folded in here: the no-bootstrap case needs it to run before
// printStartRun(), unlike the other two.
static void setIdentifiersApplyFixFactorAndPrint(DataOutputStream &ostream, StrDblMatrix &dm,
                                                  std::vector<string> &names,
                                                  bool useFixFactor, float fixfactor) {
	dm.setIdentifiers(names);
	if(useFixFactor) {
		applyFixFactor(dm,fixfactor);
	}
	ostream.print(dm);
}

// Whole-matrix output path (StrDblMatrix, fillMatrix()/print()) - the
// third of main()'s three output branches, distinct in shape from
// runRowStreamingOutput() above (different matrix type, no per-row
// streaming), so kept separate rather than forced into one function.
static void runFullMatrixOutput(DataInputStream &istream, DataOutputStream &ostream,
                                 sequence_translation_model &trans_model,
                                 int ndatasets, int numboot, bool no_incl_orig,
                                 bool useFixFactor, float fixfactor,
                                 InputFormat input_format) {
	StrDblMatrix dm;
	//open infile
	// THE DATA WE WILL PROCESS
	std::vector<Sequence> seqs;
	std::vector<string> names;
	std::vector<DNA_b128_String> b128seqs;
	Extrainfos extrainfos;

	//for each dataset in the files
	for ( int ds = 0 ; ds < ndatasets || input_format == InputFormat::Xml ; ds++ ){
		//no bootstrapping
		std::string runId;
		if ( !no_incl_orig && numboot == 0){//only need to create one distance matrix
			if ( ! istream.read(b128seqs,runId,names,extrainfos)) { break;
}
			fillMatrix(dm, b128seqs, trans_model);
			ostream.printStartRun(names,runId,extrainfos);
			//          freeXmlStrings(extrainfos);
			setIdentifiersApplyFixFactorAndPrint(ostream, dm, names, useFixFactor, fixfactor);
		}
		//bootstrapping
		else{
			//read original sequences

			if ( ! istream.readSequences(seqs,runId,extrainfos)) { break;
}

			names.clear();names.reserve(seqs.size());
			for(auto & seq : seqs) {
				names.push_back(seq.name);
}

			ostream.printStartRun(names,runId,extrainfos);
			//          freeXmlStrings(extrainfos);
			if ( !no_incl_orig ){//create the distance matrix for the original sequences
				Sequences2DNA_b128(seqs,b128seqs);
				fillMatrix(dm, b128seqs, trans_model);
				setIdentifiersApplyFixFactorAndPrint(ostream, dm, names, useFixFactor, fixfactor);
			}
			//start the bootstrapping
			//	  vector<Sequence> bootsequences;
			for ( int b = 0 ; b < numboot ; b++ ){
				//Sequence::bootstrapSequences(seqs,bootsequences);
				//Sequences2DNA_b128(bootsequences,b128seqs);
				bootstrapSequences(seqs,b128seqs);
				fillMatrix(dm, b128seqs, trans_model);
				setIdentifiersApplyFixFactorAndPrint(ostream, dm, names, useFixFactor, fixfactor);
			}
		}
		ostream.printEndRun();
	}//end data set loop

	//OUTPUT THE TREES
}////Mehmood's Changes End : email: malagori@kth.se

int
main(int argc,
		char **argv){
	setStdoutBinaryMode();

	FastdistOptions opts = parseArgs(argc, argv);
	TRY_EXCEPTION();

	handleEarlyExitFlags(opts);

	//-----------------------------------------------------
	// EVOLUTIONARY MODEL
	sequence_translation_model trans_model = buildTranslationModel(opts);

	//----------------------------------------------
	// BOOTSTRAPPING
	int numboot = opts.bootstraps;
	bool no_incl_orig = opts.no_incl_orig;

	seedRandom(opts);

	bool useFixFactor = opts.fixfactor_given;
	float fixfactor = opts.fixfactor;
	int ndatasets = opts.number_of_runs;

	//FINNISHED PARSING ARGS
	//---------------------------------------------------------
	// START BUILING MATRICES
	//

	try {
		std::unique_ptr<DataInputStream> istream = buildInputStream(opts);
		std::unique_ptr<DataOutputStream> ostream = buildOutputStream(opts);

		//Mehmood's Changes here : email: malagori@kth.se
		if (opts.output_format == OutputFormat::Binary ) {
			runRowStreamingOutput(*istream, *ostream, trans_model, ndatasets, numboot,
			                       no_incl_orig, useFixFactor, fixfactor,
			                       opts.input_format, RowOutputMode::Regular);
		} else if ( opts.memory_efficient ) {
			runRowStreamingOutput(*istream, *ostream, trans_model, ndatasets, numboot,
			                       no_incl_orig, useFixFactor, fixfactor,
			                       opts.input_format, RowOutputMode::MemoryEfficient);
		}
		else{
			runFullMatrixOutput(*istream, *ostream, trans_model, ndatasets, numboot,
			                     no_incl_orig, useFixFactor, fixfactor,
			                     opts.input_format);
		}
	}

	catch(...){
		throw;
	}

	CATCH_EXCEPTION()
	catch(...){
		cerr << "Unknown (non-Exception) error" << endl;
	}
	return 0;
}
