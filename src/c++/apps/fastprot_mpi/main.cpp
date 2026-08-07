#include "config.h"
#include "CliHelpFormatter.hpp"
#include "fastphylo/core/log_utils.hpp"
#include "fastphylo/core/file_utils.hpp"
#include "PhylipMaInputStream.hpp"
#include "DataInputStream.hpp"
#include "FastaInputStream.hpp"
#include "DataOutputStream.hpp"
#include "XmlOutputStream.hpp"
#include "fastphylo/io/fileFormatSchema.hpp"
#include "ProtDistCalc.hpp"
#include "ProtSeqUtils.hpp"
#include "fastphylo/core/DistanceMatrix.hpp"
#include "BinaryDmOutputStream.hpp"
#include "PhylipDmOutputStream.hpp"
// isatty()/STDIN_FILENO are POSIX-only (<unistd.h> doesn't exist on
// MSVC) - see the same fix in fastprot/main.cpp and fnj/main.cpp.
// fastprot_mpi is unbuilt/unverified in this environment regardless
// (no MPI installed), but fixed here too for consistency.
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
static bool stdinIsATerminal() { return _isatty(_fileno(stdin)) != 0; }
static void setStdoutBinaryMode() { _setmode(_fileno(stdout), _O_BINARY); }
#else
#include <unistd.h>
static bool stdinIsATerminal() { return isatty(STDIN_FILENO) != 0; }
static void setStdoutBinaryMode() {}
#endif

#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <numeric>

#include "mpi.h"

#ifdef WITH_LIBXML
#include "XmlInputStream.hpp"
#endif //WITH_LIBXML

// gengetopt_migration_plan.md, Phase D: fastprot_mpi is the fourth and
// last app migrated off gengetopt onto CLI11 (fnj/fastdist/fastprot
// were first, Phases B/C/D). Same FastprotOptions-shaped struct as
// fastprot's - fastprot_mpi's .ggo is a near-duplicate of fastprot's,
// minus --sd (never had it) and with a smaller --distance-function
// choice set (WAG/JTT/DAY/ARVE/MVR/LG only, no ID/JC/JCK/JCSS).
//
// UNVERIFIED: this file can't be built or run in the environment this
// migration was done in - BUILD_WITH_MPI defaults OFF and no MPI
// installation is available here, so this has never actually been
// compiled since before this change either (see the untouched
// pre-existing comment below and project_layout_plan.md's own note
// that fastprot_mpi was "never build-tested"). It also uses the
// legacy `MPI::` C++ bindings, deprecated since MPI-3 and removed
// outright from modern MPI implementations (OpenMPI 5+, recent
// MPICH) - a separate, likely-fatal, pre-existing problem this
// migration does not attempt to fix.
//
// One incidental behavior change from this migration, flagged rather
// than silently carried forward: the pre-existing code read
// args_info.input_format_arg one line *before* cmdline_parser()
// populated it (undefined-behavior read of uninitialized memory,
// harmless in practice only because WITH_LIBXML defaults ON so the
// #ifndef block doesn't normally compile in - the same bug already
// documented for fastprot/main.cpp before its own Phase D migration).
// FastprotMpiOptions's parseArgs() always fully parses before
// returning, so opts.input_format is a well-defined default
// (InputFormat::Fasta) rather than uninitialized memory if that
// #ifndef block ever does compile in - a strict improvement, not a
// behavior this migration tried to preserve bit-for-bit.
enum class InputFormat { Fasta, Phylip, Xml };
enum class OutputFormat { Phylip, Xml, Binary };

static constexpr const char *FASTPROT_MPI_VERSION = "fastprot_mpi 1.0.0";

struct FastprotMpiOptions {
	std::string outfile;
	bool outfile_given = false;
	InputFormat input_format = InputFormat::Fasta;
	bool memory_efficient = false;
	OutputFormat output_format = OutputFormat::Xml;
	int bootstraps = 0;
	bool no_incl_orig = false;
	int seed = 0;
	bool seed_given = false;
	model_type distance_function = wag;
	std::string model_file;
	bool model_file_given = false;
	bool remove_indels = false;
	bool maximum_likelihood = false;
	bool pfam = false;
	int speed = 4;
	bool print_relaxng_input = false;
	bool print_relaxng_output = false;
	std::string inputfile;
	bool inputfile_given = false;
};

static FastprotMpiOptions parseArgs(int argc, char **argv) {
	FastprotMpiOptions opts;
	CLI::App app{"Computes distance matrices out of multialignments of protein sequences", "fastprot_mpi"};
	app.formatter(std::make_shared<FastphyloHelpFormatter>());
	app.footer("\nIf FILE is not specified the input is read from stdin\n\n"
	           "Example usage of this program can be found at its home page\n"
	           "http://fastphylo.sourceforge.net/\n");
	app.set_version_flag("-V,--version", FASTPROT_MPI_VERSION);

	app.add_option("-o,--outfile", opts.outfile,
	                "output filename. If not specified, output is written to stdout")
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
	                "(possible values: phylip, xml, binary; default: xml)")
	    ->transform(CLI::CheckedTransformer(output_format_map).description(""));

	app.add_option("-b,--bootstraps", opts.bootstraps, "Bootstrap num times and create matrix for each");

	app.add_flag("-k,--no-incl-orig", opts.no_incl_orig,
	             "If the distance matrix from the original sequences should not be "
	             "included - for bootstrapping");

	CLI::Option *seed_opt = app.add_option(
	    "-R,--seed", opts.seed, "Random seed. If not specified the current timestamp will be used");

	static const std::map<std::string, model_type> distance_function_map{
	    {"WAG", wag}, {"JTT", jtt}, {"DAY", day}, {"ARVE", arve}, {"MVR", mvr}, {"LG", lg}};
	app.add_option("-D,--distance-function", opts.distance_function,
	                "Distance function (possible values: WAG, JTT, DAY, ARVE, MVR, LG; "
	                "default: WAG)")
	    ->transform(CLI::CheckedTransformer(distance_function_map).description(""));

	CLI::Option *model_file_opt =
	    app.add_option("-F,--model-file", opts.model_file,
	                    "Read matrix and equilibrium distribution from file, when used "
	                    "--distance-function is disregarded")
	        ->type_name("filename");

	app.add_flag("-i,--remove-indels", opts.remove_indels, "Remove gap columns. A gap is denoted by '-'.");

	app.add_flag("-m,--maximum-likelihood", opts.maximum_likelihood,
	             "Compute a Maximum Likelihood estimate instead. Can not be used with "
	             "--distance-function=ID, JC, JCK or JCSS");

	app.add_flag("-p,--pfam", opts.pfam,
	             "use a normal distribution as distance prior, estimated from Pfam 7.2");

	// Description text kept verbatim from the .ggo even though it says
	// "Default is 5. Valid range is [1,10]" - the .ggo's own actual
	// default="4"/values 1-8 disagree with its own description; a
	// pre-existing mismatch, not this migration's place to fix (same
	// note as fastprot's identical --speed option).
	app.add_option("-s,--speed", opts.speed,
	               "'Speed'. High speed results in low precision, only affects ED "
	               "calculations. Default is 5. Valid range is [1,10].")
	    ->check(CLI::IsMember({1, 2, 3, 4, 5, 6, 7, 8}));

	app.add_flag("-P,--print-relaxng-input", opts.print_relaxng_input,
	             "print the Relax NG schema for the XML input format (Fastphylo protein "
	             "sequence XML format) and then exit");

	app.add_flag("-w,--print-relaxng-output", opts.print_relaxng_output,
	             "print the Relax NG schema for the XML output format (Fastphylo distance "
	             "matrix XML format) and then exit.");

	CLI::Option *input_opt =
	    app.add_option("FILE", opts.inputfile, "input file. If not specified, input is read from stdin");

	try {
		app.parse(argc, argv);
	} catch (const CLI::ParseError &e) {
		app.exit(e);
		std::exit(dynamic_cast<const CLI::Success *>(&e) != nullptr ? EXIT_SUCCESS : EXIT_FAILURE);
	}

	opts.outfile_given = app["--outfile"]->count() > 0;
	opts.seed_given = seed_opt->count() > 0;
	opts.model_file_given = model_file_opt->count() > 0;
	opts.inputfile_given = input_opt->count() > 0;
	return opts;
}


int main (int argc, char **argv){
	setStdoutBinaryMode();
	if(stdinIsATerminal() && argc==1) {
	    cout<<"No input data or parameters. Use -h,--help for more information"<<endl;
	    exit(EXIT_FAILURE);
	  }

	int rank, size;
	MPI::Status status;         // status of communication

	MPI::Init(argc, argv );
	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int nr_seqs = -1;
	int seq_length = -1;
	int total_nr_dists = -1;
	int local_nr_dists = -1;

	double starttime, endtime;

	prot_sequence_translation_model trans_model;
	TRY_EXCEPTION();

	/*
  int inf = 0;
  while (inf == 0) {
    sleep(1000);
  }*/

	if (rank == 0) {
		FastprotMpiOptions opts = parseArgs(argc, argv);

#ifndef WITH_LIBXML
		if (opts.input_format == InputFormat::Xml){
			std::cerr << "The software was built with WITH_LIBXML=OFF. Please rebuild it if you want XML functionality." << std::endl; exit(EXIT_FAILURE);
		}
#endif // WITH_LIBXML

		if (opts.print_relaxng_input && opts.print_relaxng_output){
			std::cerr << "error: --print-relaxng-input and --print-relaxng-output can not be used at the same time" << std::endl; exit(EXIT_FAILURE);
		}

		if (opts.print_relaxng_input) {
			std::cout << fastphylo_prot_sequence_xml_relaxngstr << std::endl;
			exit(EXIT_SUCCESS);
		}
		if (opts.print_relaxng_output) {
			std::cout << fastphylo_distance_matrix_xml_relaxngstr << std::endl;
			exit(EXIT_SUCCESS);
		}


		//--------------------------------------------------------------
		// Read translation model

		// prot_sequence_translation_model trans_model;

		if (! opts.model_file_given){
			trans_model.model = opts.distance_function;
		} else {
			// read file from opts.model_file (never implemented, same as before this migration)
		}



		if (opts.maximum_likelihood &&
				(trans_model.model == id || trans_model.model == jc ||
						trans_model.model == jck || trans_model.model == jcss )) {
			std::cerr << "error: --maximum-likelihood can not be used with --distance-function=ID, JC, JCK or JCSS" << std::endl;
			exit(EXIT_FAILURE);
		}

		trans_model.step_size = opts.speed;

		if (opts.pfam)
			trans_model.tp = norm;
		else
			trans_model.tp = flat;

		trans_model.ml = opts.maximum_likelihood;
		bool remove_indels = opts.remove_indels;
		int ndatasets = 1;

		//----------------------------------------------
		// BOOTSTRAPPING
		int numboot = opts.bootstraps;
		bool no_incl_orig = opts.no_incl_orig;

		if ( opts.seed_given )
			srand((unsigned int )opts.seed);
		else
			srand((unsigned int)time(NULL));

		try {
			char * inputfilename = 0;
			char * outputfilename = 0;

			DataInputStream *istream;
			DataOutputStream *ostream;

			if ( opts.inputfile_given ) {
				inputfilename = const_cast<char *>(opts.inputfile.c_str());
			}

			if( opts.outfile_given )
			{  outputfilename = const_cast<char *>(opts.outfile.c_str());  }

			switch ( opts.input_format )
			{
			case InputFormat::Fasta: istream = new FastaInputStream(inputfilename);  break;
			case InputFormat::Phylip: istream = new PhylipMaInputStream(inputfilename);  break;
#ifdef WITH_LIBXML
			case InputFormat::Xml: istream = new XmlInputStream(inputfilename); break;
#endif // WITH_LIBXML
			default: exit(EXIT_FAILURE);
			}
			bool binary_format_type=opts.memory_efficient;
			switch ( opts.output_format )
			{
			case OutputFormat::Phylip: ostream = new PhylipDmOutputStream(outputfilename);  break;
			case OutputFormat::Xml: ostream = new XmlOutputStream(outputfilename); break;
			case OutputFormat::Binary:
				ostream = new BinaryDmOutputStream(outputfilename);
			        binary_format_type=true;
			        break;
			default: exit(EXIT_FAILURE);
			}

			StrDblMatrix dm;

			std::vector<Sequence> seqs;
			std::vector<std::string> names;
			Extrainfos extrainfos;

			starttime = MPI::Wtime();


			//for each dataset in the file
			for ( int ds = 0 ; ds < ndatasets || opts.input_format == InputFormat::Xml ; ds++ ){
				std::string runId("");

				if (! istream->read(seqs, runId, names, extrainfos)) break;

				// Length of null-terminated sequence
				seq_length = seqs[0].seq.size()+1;
				nr_seqs = seqs.size();

				// Send translation model to workers
				int buf[6];
				buf[0] = trans_model.model;
				buf[1] = trans_model.ml;
				buf[2] = trans_model.step_size;
				buf[3] = trans_model.tp;
				buf[4] = nr_seqs;
				buf[5] = seq_length;


				MPI::COMM_WORLD.Bcast(buf, 6, MPI::INT, 0);


				total_nr_dists = (nr_seqs*nr_seqs-nr_seqs)/2;
				local_nr_dists = total_nr_dists/size;
				if (rank < (total_nr_dists % size)) {
					local_nr_dists++;
				}

				if (remove_indels)
					remove_gaps(seqs);

				if (!no_incl_orig){
					// Copy sequences into buffer
					std::vector<char> seq_buf(seq_length*nr_seqs);
					for (int i=0; i<nr_seqs; i++)
						strcpy(&seq_buf[i*seq_length], seqs[i].seq.c_str());

					// Distribute sequences among workers
					MPI::COMM_WORLD.Bcast(&seq_buf[0], seq_length*nr_seqs, MPI::CHAR, 0);

					std::vector<double> dv(total_nr_dists);

					// Do calculations
					calculate_ed_dists_mpi(seqs, dv, trans_model, rank, size);


					// Prepare receive counts for gatherv
					std::vector<int> rec_count(size, 0);
					std::vector<int> displ(size, 0);
					int tot = 0;
					for (int i=0; i<size; i++) {
						displ[i] = tot;
						rec_count[i] = total_nr_dists / size;
						if (i < (total_nr_dists % size)) {
							rec_count[i]++;
						}
						tot += rec_count[i];
					}

					// The master process collects all the data
					MPI::COMM_WORLD.Gatherv(&dv[0], local_nr_dists, MPI::DOUBLE,
							&dv[0], &rec_count[0], &displ[0], MPI::DOUBLE, 0);


					// Copy data to distance matrix
					dm.resize(seqs.size());
					dm.setIdentifiers(names);
					int dv_ind = 0;
					for (int i=0; i<nr_seqs; i++) {
						for (int j=i+1; j<nr_seqs; j++) {
							dm.setDistance(i, j, dv[dv_ind++]);
						}
					}

					ostream->printStartRun(names, runId, extrainfos);
					ostream->printHeader(seqs.size());
					ostream->print(dm);
				}
				// Bootstrapping
				for (int b=0; b < numboot; b++){
					std::vector<Sequence> bseqs;
					bootstrap_sequences(seqs, bseqs);

					calculate_distances(bseqs, dm, trans_model);

					dm.setIdentifiers(names);
					ostream->printStartRun(names, runId, extrainfos);
					ostream->print(dm);
				}
				if (!binary_format_type){
					ostream->printEndRun();
				}
			}

			endtime = MPI::Wtime();

			delete ostream;
			delete istream;
			//printf("time: %lf\n", endtime-starttime);
		}
		catch(...){
			throw;
		}
	} else {  // Worker
		starttime = MPI::Wtime();
		// Receive translation model from master
		int buf[6];


		MPI::COMM_WORLD.Bcast(buf, 6, MPI::INT, 0);


		trans_model.model = (model_type)buf[0] ;
		trans_model.ml = buf[1];
		trans_model.step_size = buf[2];
		trans_model.tp = (type_prior)buf[3];
		nr_seqs = buf[4];
		seq_length = buf[5]; // null-terminated sequence

		total_nr_dists = (nr_seqs*nr_seqs-nr_seqs)/2;
		local_nr_dists = total_nr_dists/size;
		if (rank < (total_nr_dists % size)) {
			local_nr_dists++;
		}

		std::vector<char> seq_buf(seq_length*nr_seqs);

		// Receive sequences
		MPI::COMM_WORLD.Bcast(&seq_buf[0], seq_length*nr_seqs, MPI::CHAR, 0);

		// Copy sequences to Sequence vector
		std::vector<Sequence> seqs;
		for (int i=0; i<nr_seqs; i++) {
			seqs.push_back(Sequence("", &seq_buf[i*seq_length]));
		}

		// Calculate distances
		std::vector<double> dv(local_nr_dists);
		calculate_ed_dists_mpi(seqs, dv, trans_model, rank, size);

		// The master process collects all the data
		// receive location not used here
		MPI::COMM_WORLD.Gatherv(&dv[0], local_nr_dists, MPI::DOUBLE,
				&dv[0], 0, 0, MPI::DOUBLE, 0);
	}


	MPI::Finalize();
	CATCH_EXCEPTION();
	return 0;
}
