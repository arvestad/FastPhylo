#include "config.h"
#include "CliHelpFormatter.hpp"
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
#include "fastphylo/core/random_utils.hpp"
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <map>
// isatty()/STDIN_FILENO are POSIX-only (<unistd.h> doesn't exist on
// MSVC) - github_actions_release_builds_plan.md's Windows leg needs
// this to compile at all. _isatty()/_fileno() are MSVC's equivalents,
// declared in <io.h>.
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
static bool stdinIsATerminal()
{
    return _isatty(_fileno(stdin)) != 0;
}
// Windows' CRT defaults stdout to text mode, translating every '\n'
// this project's own code writes into "\r\n" - this project's output
// is meant to be byte-identical across platforms (RunExamples.sh byte-
// diffs it), and the code always writes '\n' explicitly, never relies
// on CRT translation, so binary mode is the correct default here, not
// just a workaround.
static void setStdoutBinaryMode()
{
    _setmode(_fileno(stdout), _O_BINARY);
}
#else
#include <unistd.h>
static bool stdinIsATerminal()
{
    return isatty(STDIN_FILENO) != 0;
}
static void setStdoutBinaryMode()
{
}
#endif

#ifdef WITH_LIBXML
#include "XmlInputStream.hpp"
#endif // WITH_LIBXML

using namespace std;

// gengetopt_migration_plan.md, Phase D: fastprot is the third app
// migrated off gengetopt onto CLI11 (fnj and fastdist were first,
// Phases B/C). FastprotOptions replaces gengetopt_args_info, same
// shape as the other two apps' Options structs. --distance-function
// binds directly to ModelMatrix.hpp's existing model_type domain enum
// via CLI::CheckedTransformer, same trick used for fnj's --method and
// fastdist's --distance-function. input-format/output-format get
// their own scoped enums (fastprot doesn't have a --number-of-runs
// option at all, unlike fastdist, so no equivalent xml-conflict check
// is needed here).
enum class InputFormat
{
    Fasta,
    Phylip,
    Xml
};
enum class OutputFormat
{
    Phylip,
    Xml,
    Binary
};

// The .ggo's own "version" line (independently pinned there, not tied
// to the top-level CMakeLists.txt's PACKAGE_VERSION="1.0.3") - kept
// exactly as-is, same reasoning as fnj's/fastdist's *_VERSION.
// github_actions_release_builds_plan.md's Phase E: PACKAGE_VERSION
// (config.h) is the single source of truth now, not an independently
// hardcoded number here.
static const std::string FASTPROT_VERSION = std::string("fastprot ") + PACKAGE_VERSION;

struct FastprotOptions
{
    std::string outfile;
    bool outfile_given = false;
    InputFormat input_format = InputFormat::Fasta;
    bool memory_efficient = false;
    OutputFormat output_format = OutputFormat::Phylip;
    int bootstraps = 0;
    bool no_incl_orig = false;
    int seed = 0;
    bool seed_given = false;
    model_type distance_function = wag;
    std::string model_file;
    bool model_file_given = false;
    bool remove_indels = false;
    bool maximum_likelihood = false;
    bool sd = false;
    bool pfam = false;
    int speed = 4;
    bool print_relaxng_input = false;
    bool print_relaxng_output = false;
    std::string inputfile;
    bool inputfile_given = false;
};

// Builds the CLI11 App and parses argv into a FastprotOptions. Exits
// the process directly on --help/--version/a parse error, same
// approach as fnj's/fastdist's parseArgs() (Phases B/C).
static FastprotOptions parseArgs(int argc, char **argv)
{
    FastprotOptions opts;
    CLI::App app{"Computes distance matrices out of multialignments of protein sequences", "fastprot"};
    app.formatter(std::make_shared<FastphyloHelpFormatter>());
    app.footer("\nIf FILE is not specified the input is read from stdin\n\n"
               "Example usage of this program can be found at its home page\n"
               "http://fastphylo.sourceforge.net/\n");
    app.set_version_flag("-V,--version", FASTPROT_VERSION);

    app.add_option("-o,--outfile", opts.outfile, "output filename. If not specified, output is written to stdout")
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

    app.add_option("-b,--bootstraps", opts.bootstraps, "Bootstrap num times and create matrix for each");

    app.add_flag("-k,--no-incl-orig", opts.no_incl_orig,
                 "If the distance matrix from the original sequences should NOT be "
                 "included - for bootstrapping");

    CLI::Option *seed_opt =
        app.add_option("-R,--seed", opts.seed, "Random seed. If not specified the current timestamp will be used");

    static const std::map<std::string, model_type> distance_function_map{{"ID", id},     {"JC", jc},   {"JCK", jck},
                                                                         {"JCSS", jcss}, {"WAG", wag}, {"JTT", jtt},
                                                                         {"DAY", day},   {"MVR", mvr}, {"LG", lg}};
    app.add_option("-D,--distance-function", opts.distance_function,
                   "Distance function (possible values: ID, JC, JCK, JCSS, WAG, JTT, DAY, "
                   "MVR, LG; default: WAG)")
        ->transform(CLI::CheckedTransformer(distance_function_map).description(""));

    CLI::Option *model_file_opt = app.add_option("-F,--model-file", opts.model_file,
                                                 "Read matrix and equilibrium distribution from file, when used "
                                                 "--distance-function is disregarded")
                                      ->type_name("filename");

    app.add_flag("-i,--remove-indels", opts.remove_indels, "Remove gap columns. A gap is denoted by '-'.");

    app.add_flag("-m,--maximum-likelihood", opts.maximum_likelihood,
                 "Compute a Maximum Likelihood estimate instead. Can not be used with "
                 "--distance-function=ID, JC, JCK or JCSS or --sd");

    app.add_flag("-S,--sd", opts.sd,
                 "Not yet implemented! Output a matrix with standard deviations after the "
                 "distance matrix. Can not be used with --distance-function=ID, JC, JCK or "
                 "JCSS or --maximum-likelihood");

    app.add_flag("-p,--pfam", opts.pfam, "use a normal distribution as distance prior, estimated from Pfam 7.2");

    // Description text kept verbatim from the .ggo even though it says
    // "Default is 5. Valid range is [1,10]" - the .ggo's own actual
    // default="4"/values 1-8 disagree with its own description; a
    // pre-existing mismatch, not this migration's place to fix.
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

    try
    {
        app.parse(argc, argv);
    }
    catch (const CLI::ParseError &e)
    {
        app.exit(e);
        // gengetopt exits 1 uniformly on any parse problem; CLI11's own
        // exit codes vary per error type. --help/--version (CLI::Success
        // and its subclasses) are the one case that's actually 0.
        std::exit(dynamic_cast<const CLI::Success *>(&e) != nullptr ? EXIT_SUCCESS : EXIT_FAILURE);
    }

    opts.outfile_given = app["--outfile"]->count() > 0;
    opts.seed_given = seed_opt->count() > 0;
    opts.model_file_given = model_file_opt->count() > 0;
    opts.inputfile_given = input_opt->count() > 0;
    return opts;
}

// Lint Phase 5 (lint_phase5_refactors.md): main() was a single 220-line
// function (cognitive complexity 59, threshold 25). Extracted below, in
// the order main() calls them.
static void handleEarlyExitFlags(const FastprotOptions &opts)
{
#ifndef WITH_LIBXML
    if (opts.input_format == InputFormat::Xml)
    {
        cerr << "The software was built with WITH_LIBXML=OFF. Please rebuild it if you want XML functionality." << endl;
        exit(EXIT_FAILURE);
    }
#endif // WITH_LIBXML
    if (opts.print_relaxng_input && opts.print_relaxng_output)
    {
        cerr << "error: --print-relaxng-input and --print-relaxng-output can not be used at the same time" << endl;
        exit(EXIT_FAILURE);
    }
    if (opts.print_relaxng_input)
    {
        cout << fastphylo_prot_sequence_xml_relaxngstr << endl;
        exit(EXIT_SUCCESS);
    }
    if (opts.print_relaxng_output)
    {
        cout << fastphylo_distance_matrix_xml_relaxngstr << endl;
        exit(EXIT_SUCCESS);
    }
}

static prot_sequence_translation_model buildTranslationModel(const FastprotOptions &opts)
{
    prot_sequence_translation_model trans_model;
    if (!opts.model_file_given)
    {
        trans_model.model = opts.distance_function;
    }
    else
    {
        // read file from opts.model_file (never implemented, same as before this migration)
    }
    if (opts.maximum_likelihood && (trans_model.model == id || trans_model.model == jc || trans_model.model == jck ||
                                    trans_model.model == jcss || opts.sd))
    {
        cerr << "error: --maximum-likelihood can not be used with --distance-function=ID, JC, JCK or JCSS or --sd"
             << endl;
        exit(EXIT_FAILURE);
    }
    if (opts.sd && (trans_model.model == id || trans_model.model == jc || trans_model.model == jck ||
                    trans_model.model == jcss || opts.maximum_likelihood))
    {
        cerr << "error: --sd can not be used with --distance-function=ID, JC, JCK or JCSS or --maximum-likelihood"
             << endl;
        exit(EXIT_FAILURE);
    }
    if (opts.sd)
    {
        // Remove this two lines when sd is working again
        cerr << "Aborting: The std dev feature is not yet implemented." << endl;
        exit(1);
    }
    trans_model.step_size = opts.speed;
    if (opts.pfam)
    {
        // Explicit ::norm (not std::norm, ambiguous here since this
        // file's `using namespace std;` plus Eigen pulling in
        // <complex> transitively via Matrix.hpp makes the unqualified
        // name ambiguous) - this is ExpectedDistance.hpp's type_prior
        // enum value, unrelated to std::norm's complex-number meaning.
        trans_model.tp = ::norm;
    }
    else
    {
        trans_model.tp = flat;
    }
    trans_model.sd = opts.sd;
    trans_model.ml = opts.maximum_likelihood;
    return trans_model;
}

static std::unique_ptr<DataInputStream> buildInputStream(const FastprotOptions &opts)
{
    char *inputfilename = nullptr;
    if (opts.inputfile_given)
    {
        inputfilename = const_cast<char *>(opts.inputfile.c_str());
    }
    switch (opts.input_format)
    {
    case InputFormat::Fasta:
        return std::make_unique<FastaInputStream>(inputfilename);
    case InputFormat::Phylip:
        return std::make_unique<PhylipMaInputStream>(inputfilename);
#ifdef WITH_LIBXML
    case InputFormat::Xml:
        return std::make_unique<XmlInputStream>(inputfilename);
#endif // WITH_LIBXML
    default:
        exit(EXIT_FAILURE);
    }
}

static std::unique_ptr<DataOutputStream> buildOutputStream(const FastprotOptions &opts)
{
    char *outputfilename = nullptr;
    if (opts.outfile_given)
    {
        outputfilename = const_cast<char *>(opts.outfile.c_str());
    }
    switch (opts.output_format)
    {
    case OutputFormat::Phylip:
        return std::make_unique<PhylipDmOutputStream>(outputfilename);
    case OutputFormat::Xml:
        return std::make_unique<XmlOutputStream>(outputfilename);
    case OutputFormat::Binary:
        // binary_dm_format_plan.md: matricesPerRun is a single value
        // fixed for the whole invocation here - ndatasets is always 1
        // for non-XML input (see processRuns() below), and numboot/
        // no_incl_orig never change per dataset even for XML input.
        return std::make_unique<BinaryDmOutputStream>(outputfilename, (opts.no_incl_orig ? 0 : 1) +
                                                                          static_cast<size_t>(opts.bootstraps));
    default:
        exit(EXIT_FAILURE);
    }
}

// Shared by processRuns()'s original-data and bootstrap-replicate
// cases below: pick ml/sd/non-sd calculate_distances(), then stamp
// identifiers - everything else around it (whether printStartRun()/
// printHeader() run first) differs between the two call sites, so
// only this common part is factored out.
//
// ml_decomp is built once, outside processRuns()'s dataset/bootstrap
// loop (see main() below) - the ML branch here calls
// calculate_ml_dists() directly with it rather than routing through
// calculate_distances(), which would rebuild the (now much more
// expensive to skip) decomposition on every single bootstrap
// replicate - see build_ml_decomposition()'s doc comment
// (ProtDistCalc.hpp) for why that matters. Non-null iff
// trans_model.ml is true (main() guarantees this).
static void calculateDistances(prot_sequence_translation_model &trans_model, vector<Sequence> &seqs, StrDblMatrix &dm,
                               StrDblMatrix &sdm, vector<string> &names, const MatrixExpm *ml_decomp)
{
    if (trans_model.ml)
    {
        calculate_ml_dists(seqs, dm, *ml_decomp);
    }
    else if (trans_model.sd)
    {
        calculate_distances(seqs, dm, trans_model, sdm);
    }
    else
    {
        calculate_distances(seqs, dm, trans_model);
    }
    dm.setIdentifiers(names);
}

static void printDistances(DataOutputStream &ostream, prot_sequence_translation_model &trans_model, StrDblMatrix &dm,
                           StrDblMatrix &sdm, bool binary_format_type)
{
    ostream.print(dm);
    if (trans_model.sd && !binary_format_type)
    {
        ostream.printSD(sdm);
    }
}

// The dataset/bootstrap run loop - reads each dataset, computes and
// prints its distance matrix, then numboot bootstrap replicates of it.
// ml_decomp: see calculateDistances()'s doc comment above.
static void processRuns(DataInputStream &istream, DataOutputStream &ostream,
                        prot_sequence_translation_model &trans_model, int ndatasets, int numboot, bool no_incl_orig,
                        bool remove_indels, bool binary_format_type, InputFormat input_format,
                        const MatrixExpm *ml_decomp)
{
    StrDblMatrix dm;
    StrDblMatrix sdm;
    vector<Sequence> seqs;
    vector<string> names;
    Extrainfos extrainfos;

    for (int ds = 0; ds < ndatasets || input_format == InputFormat::Xml; ds++)
    {
        string runId;
        if (!istream.read(seqs, runId, names, extrainfos))
        {
            break;
        }
        if (remove_indels)
        {
            remove_gaps(seqs);
        }
        if (!no_incl_orig)
        {
            calculateDistances(trans_model, seqs, dm, sdm, names, ml_decomp);
            ostream.printStartRun(names, runId, extrainfos);
            ostream.printHeader(seqs.size());
            printDistances(ostream, trans_model, dm, sdm, binary_format_type);
        }
        // Bootstrapping
        for (int b = 0; b < numboot; b++)
        {
            vector<Sequence> bseqs;
            bootstrap_sequences(seqs, bseqs);
            calculateDistances(trans_model, bseqs, dm, sdm, names, ml_decomp);
            printDistances(ostream, trans_model, dm, sdm, binary_format_type);
        }
        if (!binary_format_type)
        {
            ostream.printEndRun();
        }
    }
}

int main(int argc, char **argv)
{
    setStdoutBinaryMode();
    if (stdinIsATerminal() && argc == 1)
    {
        cout << "No input data or parameters. Use -h,--help for more information" << endl;
        exit(EXIT_FAILURE);
    }
    FastprotOptions opts = parseArgs(argc, argv);
    try
    {
        handleEarlyExitFlags(opts);

        prot_sequence_translation_model trans_model = buildTranslationModel(opts);
        bool remove_indels = opts.remove_indels;
        int ndatasets = 1;
        //----------------------------------------------
        // BOOTSTRAPPING
        int numboot = opts.bootstraps;
        bool no_incl_orig = opts.no_incl_orig;
        if (numboot == 0 && no_incl_orig)
        {
            cerr << "error: No output. No bootstrap or original data will be written" << endl;
            exit(EXIT_FAILURE);
        }
        if (opts.seed_given)
        {
            seed_rng(static_cast<unsigned int>(opts.seed));
        }
        else
        {
            // Not a security context - bootstrap resampling just needs
            // an arbitrary starting point when the user hasn't asked
            // for a reproducible one via --seed.
            seed_rng(static_cast<unsigned int>(time(nullptr)));
        }
        std::unique_ptr<DataInputStream> istream = buildInputStream(opts);
        std::unique_ptr<DataOutputStream> ostream = buildOutputStream(opts);
        bool binary_format_type = opts.memory_efficient || (opts.output_format == OutputFormat::Binary);

        // Built once, before any dataset/bootstrap replicate is
        // processed - trans_model.model never changes over the course
        // of a run, so this is the one point that's genuinely "before
        // any distance estimation happens" for the whole program, not
        // just for one dataset - see calculateDistances()'s doc
        // comment above for why this matters.
        std::optional<MatrixExpm> ml_decomp;
        if (trans_model.ml)
        {
            ml_decomp = build_ml_decomposition(trans_model.model);
        }

        processRuns(*istream, *ostream, trans_model, ndatasets, numboot, no_incl_orig, remove_indels,
                    binary_format_type, opts.input_format, ml_decomp ? &*ml_decomp : nullptr);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
    catch (...)
    {
        std::cerr << "Unknown (non-Exception) error" << std::endl;
        exit(EXIT_FAILURE);
    }
    return 0;
}
