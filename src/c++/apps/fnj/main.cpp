#include <memory>
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <unistd.h>

#include "config.h"
#include "fastphylo/core/file_utils.hpp"
#include "fastphylo/dna/NeighborJoining.hpp"
#include <ctime>
#include <ctime>
#include "fastphylo/core/log_utils.hpp"
#include "fnj_gengetopt.h"

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

int main (int argc, char **argv) {
    if((isatty(STDIN_FILENO) != 0) && argc==1) {
      cout<<"No input data or parameters. Use -h,--help for more information"<<endl;
      exit(EXIT_FAILURE);
    }
	gengetopt_args_info args_info;
	TRY_EXCEPTION();
	if (cmdline_parser (argc, argv, &args_info) != 0) {
		exit(EXIT_FAILURE);
}
#ifndef WITH_LIBXML
	if ( args_info.input_format_arg == input_format_arg_xml ) {
		cerr << "The software was built with WITH_LIBXML=OFF. Please rebuild it if you want XML functionality." << endl;
		exit(EXIT_FAILURE);
	}
#endif // WITH_LIBXML
	if ( (args_info.print_relaxng_input_given != 0u) && (args_info.print_relaxng_output_given != 0u) ) {
		cerr << "error: --print-relaxng-input and --print-relaxng-output can not be used at the same time" << endl;
		exit(EXIT_FAILURE);
	}
	if ( args_info.print_relaxng_input_given != 0u ) {
		cout << fastphylo_distance_matrix_xml_relaxngstr << std::endl;
		exit(EXIT_SUCCESS);
	};
	if ( args_info.print_relaxng_output_given != 0u ) {
		cout << fastphylo_tree_count_xml_relaxngstr << std::endl;
		exit(EXIT_SUCCESS);
	};
	//----------------------------------------------
	// DISTANCE METHODS
	std::vector<NJ_method> methods;
	if( (args_info.number_of_runs_given != 0u) && args_info.input_format_arg == input_format_arg_xml ) {
		cerr << "error: --number-of-runs can not be used together with input format xml." << endl;
		exit(EXIT_FAILURE);
	}
	switch ( args_info.method_arg ) {
	case method_arg_NJ:
		methods.push_back(NJ);
		break;
	case method_arg_FNJ:
		methods.push_back(FNJ);
		break;
	case method_arg_BIONJ:
		methods.push_back(BIONJ);
		break;
	default:
		cerr << "error: method chosen not available" << endl;
		exit(EXIT_FAILURE);
	}
	bool printCounts = args_info.print_counts_flag != 0;
	try {
		char * inputfilename = nullptr;
		char * outputfilename = nullptr;
		std::unique_ptr<DataInputStream> istream;
		std::unique_ptr<DataOutputStream> ostream;

		switch( args_info.inputs_num ) {
			case 0:
				break; /* inputfilename will be null and indicate stdin as input */
		case 1:
			inputfilename =  args_info.inputs[0];
			break;
		default: cerr << "Error: you can at most specify one input filename" << endl;
			exit(EXIT_FAILURE);
		}
		if( args_info.outfile_given != 0u ) {
			outputfilename = args_info.outfile_arg;
}
		switch ( args_info.input_format_arg ) {
			case input_format_arg_phylip:
				istream = std::make_unique<PhylipDmInputStream>(inputfilename);
				break;
			case input_format_arg_binary: istream = std::make_unique<BinaryInputStream>(inputfilename);
				break;
#ifdef WITH_LIBXML
			case input_format_arg_xml: istream = std::make_unique<XmlInputStream>(inputfilename);
				break;
#endif // WITH_LIBXML
			default:
				exit(EXIT_FAILURE);
		}
		switch (args_info.output_format_arg) {
			case output_format_arg_newick:
				ostream = std::make_unique<TreeTextOutputStream>(outputfilename);
				break;
			case output_format_arg_xml:
				ostream = std::make_unique<XmlOutputStream>(outputfilename);
				break;
			default:
				exit(EXIT_FAILURE);
		}
		//printf("%d\n", args_info.input_format_arg);
		// THE DATA WE WILL PROCESS
		vector<Sequence> seqs;
		vector<std::string> names;
		vector<DNA_b128_String> b128seqs;
		Extrainfos extrainfos;
		bool latestReadSuccessful = true;
		vector<string> speciesnames;
		readstatus status;
		int run = 0;
		status = END_OF_RUN;

		while (status == END_OF_RUN && (args_info.input_format_arg == input_format_arg_xml || run<args_info.number_of_runs_arg)) {
			string runId;
			run++;
			tree2int_map tree2count(static_cast<size_t>(args_info.bootstraps_arg * 1.3));
			str2int_hashmap name2id;
			if (args_info.input_format_arg==input_format_arg_binary) {
				StrDblMatrix dm;
				for (int runNo=1; (status = istream->readDM(dm, names, runId, extrainfos))==DM_READ; runNo++) {
					if (args_info.analyze_run_number_given != 0u) {
						if (runNo<args_info.analyze_run_number_arg) {
							continue;
}
						if (runNo>args_info.analyze_run_number_arg) {
							status=END_OF_RUN;
							break;
							}
						}
					for(size_t namei=0; namei<dm.getSize(); namei++) {
						name2id[dm.getIdentifier(namei)] = namei;
}
					buildTrees(dm, tree2count, methods,name2id);
				}
			}
			else {
				StrDblMatrix dm;
				for (int runNo=1; (status = istream->readDM(dm, names, runId, extrainfos))==DM_READ; runNo++) {
					if (args_info.analyze_run_number_given != 0u) {
						if (runNo<args_info.analyze_run_number_arg) {
							continue;
}
						if (runNo>args_info.analyze_run_number_arg) {
							status=END_OF_RUN;
							break;
							}
						}
					for(size_t namei=0; namei<dm.getSize(); namei++) {
					     name2id[dm.getIdentifier(namei)] = namei;
					}
					buildTrees(dm, tree2count, methods,name2id);
				}
			}
			if (status==END_OF_RUN) {
				ostream->print(tree2count,printCounts, runId, names, extrainfos);
}
			if (args_info.analyze_run_number_given != 0u) {
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
	cmdline_parser_free(&args_info);
	return 0;
}
