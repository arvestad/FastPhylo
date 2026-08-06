#include <memory>
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <map>
#include <unistd.h>

#include "config.h"
#include "fastphylo/core/file_utils.hpp"
#include "fastphylo/dna/NeighborJoining.hpp"
#include <ctime>
#include "fastphylo/core/log_utils.hpp"
#include "CliHelpFormatter.hpp"

#include "DataInputStream.hpp"
#include "DataOutputStream.hpp"
#include "TreeTextOutputStream.hpp"
#include "XmlOutputStream.hpp"
#include "fastphylo/io/Extrainfos.hpp"
#include "fastphylo/io/fileFormatSchema.hpp"
#include "PhylipDmInputStream.hpp"
#include "BinaryInputStream.hpp"

#ifdef WITH_LIBXML
#include "XmlInputStream.hpp"
#endif // WITH_LIBXML

using namespace std;

// gengetopt_migration_plan.md, Phase B: fnj is the first app migrated
// off gengetopt onto CLI11. FnjOptions replaces gengetopt_args_info as
// the type every helper function below takes - same shape (a plain
// struct of parsed values plus "_given" flags for the handful of
// options whose *presence*, not just their value, changes behavior).
// input-format/output-format get their own scoped enums; --method
// reuses NeighborJoining.hpp's existing NJ_method domain enum directly
// rather than inventing a parallel one, since CLI::CheckedTransformer
// can parse straight into it.
enum class InputFormat { Phylip, Xml, Binary };
enum class OutputFormat { Newick, Xml };

// The .ggo's own "version" line (independently pinned there, not tied
// to the top-level CMakeLists.txt's PACKAGE_VERSION="1.0.3") - kept
// exactly as-is rather than "fixed" to match PACKAGE_VERSION, since
// that would be a behavior change outside this migration's scope.
static constexpr const char *FNJ_VERSION = "fnj 1.0.10";

struct FnjOptions {
	std::string outfile;
	bool outfile_given = false;
	InputFormat input_format = InputFormat::Xml;
	OutputFormat output_format = OutputFormat::Xml;
	bool print_counts = false;
	int analyze_run_number = 0;
	bool analyze_run_number_given = false;
	NJ_method method = FNJ;
	int number_of_runs = 1;
	bool number_of_runs_given = false;
	int bootstraps = 0;
	bool print_relaxng_input = false;
	bool print_relaxng_output = false;
	std::string inputfile;
	bool inputfile_given = false;
};

// Builds the CLI11 App and parses argv into a FnjOptions. Exits the
// process directly on --help/--version/a parse error, same as
// gengetopt's cmdline_parser() did - kept as its own function so
// main() reads the same as it did pre-migration.
static FnjOptions parseArgs(int argc, char **argv) {
	FnjOptions opts;
	CLI::App app{"Builds phylogenetic trees", "fnj"};
	app.formatter(std::make_shared<FastphyloHelpFormatter>());
	app.footer("\nExample usage of this program can be found at its home page\n"
	           "http://fastphylo.sourceforge.net/\n");
	app.set_version_flag("-V,--version", FNJ_VERSION);

	app.add_option("-o,--outfile", opts.outfile,
	                "output filename. If not specifed, output is written to stdout")
	    ->type_name("filename");

	static const std::map<std::string, InputFormat> input_format_map{
	    {"phylip", InputFormat::Phylip}, {"xml", InputFormat::Xml}, {"binary", InputFormat::Binary}};
	app.add_option("-I,--input-format", opts.input_format,
	                "input format. 'xml' means the 'Fastphylo distance matrix XML format' "
	                "(possible values: phylip, xml, binary; default: xml)")
	    ->transform(CLI::CheckedTransformer(input_format_map).description(""));

	static const std::map<std::string, OutputFormat> output_format_map{
	    {"newick", OutputFormat::Newick}, {"xml", OutputFormat::Xml}};
	app.add_option("-O,--output-format", opts.output_format,
	                "output format. 'xml' means the 'Fastphylo tree count XML format' "
	                "(possible values: newick, xml; default: xml)")
	    ->transform(CLI::CheckedTransformer(output_format_map).description(""));

	app.add_flag("-c,--print-counts", opts.print_counts,
	             "print the tree count before each the newick tree. This flag has no "
	             "effect on the XML output format.");

	app.add_option("-a,--analyze-run-number", opts.analyze_run_number,
	               "Determines which dataset should be analyzed with 1 being the first "
	               "dataset. By default all are analyzed");

	static const std::map<std::string, NJ_method> method_map{
	    {"NJ", NJ}, {"FNJ", FNJ}, {"BIONJ", BIONJ}};
	app.add_option("-m,--method", opts.method,
	               "reconstruction method to apply (possible values: NJ, FNJ, BIONJ; default: FNJ)")
	    ->transform(CLI::CheckedTransformer(method_map).description(""));

	// Accepted for CLI compatibility only - never consumed by fnj's own
	// code, confirmed by grep before this migration; not this
	// migration's place to silently drop or start honoring it.
	int dm_per_run = 1;
	app.add_option("-d,--dm-per-run", dm_per_run,
	               "nr of Distance matrices per run. Is only used if the input format is phylip");

	app.add_option("-r,--number-of-runs", opts.number_of_runs,
	               "nr of runs. Is only used if the input format is phylip");

	app.add_option("-b,--bootstraps", opts.bootstraps, "number of boot straps");

	app.add_flag("-p,--print-relaxng-input", opts.print_relaxng_input,
	             "print the Relax NG schema for the XML input format (Fastphylo "
	             "distance matrix XML format) and then exit");

	app.add_flag("-w,--print-relaxng-output", opts.print_relaxng_output,
	             "print the Relax NG schema for the XML output format (Fastphylo tree "
	             "count XML format) and then exit.");

	CLI::Option *input_opt =
	    app.add_option("FILE", opts.inputfile, "input file. If not specified, input is read from stdin");

	try {
		app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
		app.exit(e);
		// gengetopt exits 1 uniformly on any parse problem; CLI11's own
		// exit codes vary per error type (e.g. 105/109). --help/--version
		// (CLI::Success and its subclasses) are the one case that's
		// actually 0, matching cmdline_parser()'s own -h/-V handling.
		std::exit(dynamic_cast<const CLI::Success *>(&e) != nullptr ? EXIT_SUCCESS : EXIT_FAILURE);
	}

	opts.outfile_given = app["--outfile"]->count() > 0;
	opts.analyze_run_number_given = app["--analyze-run-number"]->count() > 0;
	opts.number_of_runs_given = app["--number-of-runs"]->count() > 0;
	opts.inputfile_given = input_opt->count() > 0;
	return opts;
}

//
// Builds trees using the supplied distance methods.
// Each created tree is added to the tree2count map.
//

template<class T> void buildTrees(T &dm, tree2int_map &tree2count, std::vector<NJ_method> &methods, str2int_hashmap &name2id) {
	SequenceTree tree;

	for(auto & method : methods){
		computeNJTree(dm,tree,method);
		tree.makeCanonical(name2id);
		auto iter = tree2count.find(tree);
		if(iter!=tree2count.end()) {
			iter->second++;
		} else {
			tree2count[tree] = 1;
}
	}
}

// Lint Phase 5 (lint_phase5_refactors.md): main() was a single 162-line
// function (cognitive complexity 69, threshold 25). handleEarlyExitFlags()
// covers every check below except the isatty() one (which has to run
// before argv is parsed at all, so it stays inline in main()).
static void handleEarlyExitFlags(const FnjOptions &opts) {
#ifndef WITH_LIBXML
	if ( opts.input_format == InputFormat::Xml ) {
		cerr << "The software was built with WITH_LIBXML=OFF. Please rebuild it if you want XML functionality." << endl;
		exit(EXIT_FAILURE);
	}
#endif // WITH_LIBXML
	if ( opts.print_relaxng_input && opts.print_relaxng_output ) {
		cerr << "error: --print-relaxng-input and --print-relaxng-output can not be used at the same time" << endl;
		exit(EXIT_FAILURE);
	}
	if ( opts.print_relaxng_input ) {
		cout << fastphylo_distance_matrix_xml_relaxngstr << std::endl;
		exit(EXIT_SUCCESS);
	}
	if ( opts.print_relaxng_output ) {
		cout << fastphylo_tree_count_xml_relaxngstr << std::endl;
		exit(EXIT_SUCCESS);
	}
	if( opts.number_of_runs_given && opts.input_format == InputFormat::Xml ) {
		cerr << "error: --number-of-runs can not be used together with input format xml." << endl;
		exit(EXIT_FAILURE);
	}
}

static std::unique_ptr<DataInputStream> buildInputStream(const FnjOptions &opts) {
	char *inputfilename = nullptr;
	if ( opts.inputfile_given ) {
		inputfilename = const_cast<char *>(opts.inputfile.c_str());
	}
	switch ( opts.input_format ) {
		case InputFormat::Phylip:
			return std::make_unique<PhylipDmInputStream>(inputfilename);
		case InputFormat::Binary:
			return std::make_unique<BinaryInputStream>(inputfilename);
#ifdef WITH_LIBXML
		case InputFormat::Xml:
			return std::make_unique<XmlInputStream>(inputfilename);
#endif // WITH_LIBXML
		default:
			exit(EXIT_FAILURE);
	}
}

static std::unique_ptr<DataOutputStream> buildOutputStream(const FnjOptions &opts) {
	char *outputfilename = nullptr;
	if( opts.outfile_given ) {
		outputfilename = const_cast<char *>(opts.outfile.c_str());
	}
	switch (opts.output_format) {
		case OutputFormat::Newick:
			return std::make_unique<TreeTextOutputStream>(outputfilename);
		case OutputFormat::Xml:
			return std::make_unique<XmlOutputStream>(outputfilename);
		default:
			exit(EXIT_FAILURE);
	}
}

// One run's worth of distance matrices, read and turned into trees.
// input_format_arg_binary used to get its own hand-duplicated copy of
// this loop, byte-for-byte identical to the non-binary case below it -
// the input format only matters to buildInputStream()'s choice of
// DataInputStream implementation above, never to this loop itself.
static void processRuns(DataInputStream &istream, tree2int_map &tree2count,
                         std::vector<NJ_method> &methods, str2int_hashmap &name2id,
                         std::vector<std::string> &names, std::string &runId, Extrainfos &extrainfos,
                         const FnjOptions &opts, readstatus &status) {
	StrDblMatrix dm;
	for (int runNo=1; (status = istream.readDM(dm, names, runId, extrainfos))==DM_READ; runNo++) {
		if (opts.analyze_run_number_given) {
			if (runNo<opts.analyze_run_number) {
				continue;
			}
			if (runNo>opts.analyze_run_number) {
				status=END_OF_RUN;
				break;
			}
		}
		for(size_t namei=0; namei<dm.getSize(); namei++) {
			name2id[dm.getIdentifier(namei)] = namei;
		}
		buildTrees(dm, tree2count, methods, name2id);
	}
}

int main (int argc, char **argv) {
    if((isatty(STDIN_FILENO) != 0) && argc==1) {
      cout<<"No input data or parameters. Use -h,--help for more information"<<endl;
      exit(EXIT_FAILURE);
    }
	FnjOptions opts = parseArgs(argc, argv);
	TRY_EXCEPTION();
	handleEarlyExitFlags(opts);
	std::vector<NJ_method> methods{opts.method};
	bool printCounts = opts.print_counts;
	try {
		std::unique_ptr<DataInputStream> istream = buildInputStream(opts);
		std::unique_ptr<DataOutputStream> ostream = buildOutputStream(opts);

		vector<std::string> names;
		Extrainfos extrainfos;
		readstatus status = END_OF_RUN;
		int run = 0;

		while (status == END_OF_RUN && (opts.input_format == InputFormat::Xml || run<opts.number_of_runs)) {
			string runId;
			run++;
			tree2int_map tree2count(static_cast<size_t>(opts.bootstraps * 1.3));
			str2int_hashmap name2id;

			processRuns(*istream, tree2count, methods, name2id, names, runId, extrainfos, opts, status);

			if (status==END_OF_RUN) {
				ostream->print(tree2count,printCounts, runId, names, extrainfos);
}
			if (opts.analyze_run_number_given) {
				break;
}
		}//end run loop
	}
	catch(...){
		throw;
	}
	CATCH_EXCEPTION()
	catch(...){
		std::cerr << "Unknown (non-Exception) error" << std::endl;
	}
	return 0;
}
