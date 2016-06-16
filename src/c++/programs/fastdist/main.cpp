///////////////////////////////////////////////
//                                           //
// File: fastdist.cpp                        //
//                                           //
// Created on: <14-Mar-2005 14:21:14 cherub> //
//                                           //
// Author: Isaac Elias,Mehmood Alam Khan     //
// Email: isaac@nada.kth.se, malagori@kth.se //
///////////////////////////////////////////////

#include "Sequences2DistanceMatrix.hpp"

#include <string>
#include <iostream>
#include <time.h>
#include <fstream>
#include <assert.h>

#include "config.h"
#include "file_utils.hpp"
#include <iomanip>
#include "log_utils.hpp"
#include "BinaryDmOutputStream.hpp"
#include "fastdist_gengetopt.h"
#include "NeighborJoining.hpp"
#include "DataInputStream.hpp"
#include "PhylipMaInputStream.hpp"
#include "FastaInputStream.hpp"
#include "DataOutputStream.hpp"
#include "Extrainfos.hpp"
#include "fileFormatSchema.hpp"
#include "XmlOutputStream.hpp"
#include "PhylipDmOutputStream.hpp"

#ifdef WITH_LIBXML
#include "XmlInputStream.hpp"
#endif // WITH_LIBXML

using namespace std;

int
main(int argc,
		char **argv){

	gengetopt_args_info args_info;
	TRY_EXCEPTION();

	sequence_translation_model trans_model;



#ifndef WITH_LIBXML
	if ( args_info.input_format_arg == input_format_arg_xml ) {
		cerr << "The software was built with WITH_LIBXML=OFF. Please rebuild it if you want XML functionality." << endl; exit(EXIT_FAILURE);
	}
#endif // WITH_LIBXML

	if (cmdline_parser (argc, argv, &args_info) != 0)
		exit(EXIT_FAILURE);

	if ( args_info.print_relaxng_input_given && args_info.print_relaxng_output_given ) {
		cerr << "error: --print-relaxng-input and --print-relaxng-output can not be used at the same time" << endl; exit(EXIT_FAILURE);
	}

	if ( args_info.print_relaxng_input_given ) {  cout << fastphylo_sequence_xml_relaxngstr << std::endl;  exit(EXIT_SUCCESS);   };
	if ( args_info.print_relaxng_output_given ) {  cout << fastphylo_distance_matrix_xml_relaxngstr << std::endl;  exit(EXIT_SUCCESS);   };

	if ( args_info.number_of_runs_given && args_info.input_format_arg != input_format_arg_phylip  ) {
		cerr << "error: --number-of-runs can only be used together with --input-format=phylip " << endl; exit(EXIT_FAILURE);
	}

	//-----------------------------------------------------
	// EVOLUTIONARY MODEL

	switch ( args_info.distance_function_arg )
	{
	case distance_function_arg_JC : trans_model.model = JC; break;
	case distance_function_arg_K2P : trans_model.model = K2P; break;
	case distance_function_arg_TN93 : trans_model.model = TN93; break;
	case distance_function_arg_HAMMING : trans_model.model = HAMMING_DISTANCE; break;
	default: cerr << "error: model chosen not available" << endl; exit(EXIT_FAILURE);
	}

	trans_model.no_tstvratio = args_info.no_tstvratio_given;
	trans_model.tstvratio = args_info.tstvratio_arg;
	trans_model.pyrtvratio =  args_info.pyrtvratio_arg;

	//----------------------------------------------
	// BOOTSTRAPPING
	int numboot = args_info.bootstraps_arg;
	bool no_incl_orig = args_info.no_incl_orig_given;

	if ( args_info.seed_given )
		srand((unsigned int )args_info.seed_arg);
	else
		srand((unsigned int)time(NULL));

	//-----------------------------------------------
	// AMBIGUITIES

	trans_model.no_ambiguities = args_info.no_ambiguities_given;
	trans_model.no_ambig_resolve = args_info.no_ambig_resolve_given;
	trans_model.no_transition_probs = args_info.no_transprob_given;

	switch ( args_info.ambiguity_frequency_model_arg )
	{
	case ambiguity_frequency_model_arg_UNI : trans_model.use_base_freqs = false; break;
	case ambiguity_frequency_model_arg_BASE : trans_model.use_base_freqs = true; break;
	default: cerr << "programming error 2..." << endl; exit(EXIT_FAILURE);
	}
	bool useFixFactor = args_info.fixfactor_given;
	float fixfactor=args_info.fixfactor_arg;
	int ndatasets = args_info.number_of_runs_arg;

	//FINNISHED PARSING ARGS
	//---------------------------------------------------------
	// START BUILING MATRICES
	//

	try {
		char * inputfilename = NULL;
		char * outputfilename = NULL;

		DataInputStream *istream;
		DataOutputStream *ostream;

		switch( args_info.inputs_num )
		  {
		  case 0: break; /* inputfilename will be null and indicate stdin as input */
		  case 1: inputfilename =  args_info.inputs[0]; break;
		  default: cerr << "Error: you can at most specify one input filename" << endl; exit(EXIT_FAILURE);
		}

		if( args_info.outfile_given )
		{  outputfilename = args_info.outfile_arg;  }

		switch ( args_info.input_format_arg )
		{
		case input_format_arg_fasta: istream = new FastaInputStream(inputfilename);  break;
		case input_format_arg_phylip : istream = new PhylipMaInputStream(inputfilename);  break;
#ifdef WITH_LIBXML
		case input_format_arg_xml: istream = new XmlInputStream(inputfilename); break;
#endif // WITH_LIBXML
		default: exit(EXIT_FAILURE);
		}

		switch ( args_info.output_format_arg )
		{
		case output_format_arg_phylip: ostream = new PhylipDmOutputStream(outputfilename);  break;
		case output_format_arg_xml: ostream = new XmlOutputStream(outputfilename); break;
		//Mehmood's Changes here : email: malagori@kth.se
		case output_format_arg_binary: ostream = new BinaryDmOutputStream(outputfilename); break;
		default: exit(EXIT_FAILURE);
		}

		//Mehmood's Changes here : email: malagori@kth.se
		if (args_info.output_format_arg == output_format_arg_binary ) {
			StrFloRow dm;
			//open infile

			// THE DATA WE WILL PROCESS
			std::vector<Sequence> seqs;
			std::vector<string> names;
			std::vector<DNA_b128_String> b128seqs;

			Extrainfos extrainfos;

			//for each dataset in the files


			for ( int ds = 0 ; ds < ndatasets || args_info.input_format_arg == input_format_arg_xml ; ds++ ){
				//no bootstrapping
				std::string runId("");
				if ( !no_incl_orig && numboot == 0){//only need to create one distance matrix
					if ( ! istream->read(b128seqs,runId,names,extrainfos)) {
						break;
					}

					const size_t numberOfSequences = b128seqs.size();
					ostream->printStartRun(names,runId,extrainfos);
					ostream->printHeader(numberOfSequences);

					for(size_t i = 0; i < numberOfSequences; ++i){
						fillMatrixRow(dm, b128seqs, trans_model, i, false);
						dm.setIdentifier(names.at(i));
						if(useFixFactor) applyFixFactorRow(dm,fixfactor);
						ostream->printRow(dm, names.at(i), i, false);
					}
				}
				//bootstrapping
				else{
					//TODO: implement bootstraping with the improved version of FastDist
					StrFloRow dm;
					//read original sequences
					if ( ! istream->readSequences(seqs,runId,extrainfos)) break;
					names.clear();names.reserve(seqs.size());
					for( size_t i=0;i<seqs.size();i++)
					{
						names.push_back(seqs[i].name);
					}
					const size_t numberOfSequences = seqs.size();
					ostream->printStartRun(names,runId,extrainfos);
					ostream->printHeader(numberOfSequences);
					if ( !no_incl_orig ){//create the distance matrix for the original sequences
						Sequences2DNA_b128(seqs,b128seqs);
						for(size_t i = 0; i < numberOfSequences; ++i){
							fillMatrixRow(dm, b128seqs, trans_model, i, false);
							dm.setIdentifier(names.at(i));
							if(useFixFactor) applyFixFactorRow(dm,fixfactor);
							ostream->printRow(dm, names.at(i), i, false);
						}
					}
					//start the bootstrapping
					//	  vector<Sequence> bootsequences;
					for ( int b = 0 ; b < numboot ; b++ ){
						bootstrapSequences(seqs,b128seqs);
						std::cout << numberOfSequences << std::endl;
						ostream->printBootstrapSpliter(numberOfSequences);
						for(size_t i = 0; i < numberOfSequences; ++i){
							fillMatrixRow(dm, b128seqs, trans_model, i, false);
							dm.setIdentifier(names.at(i));
							if(useFixFactor) applyFixFactorRow(dm,fixfactor);
							ostream->printRow(dm, names.at(i), i,false);
						}
					}
				}
				ostream->printEndRun();
			}//end data set loop
		} else if ( args_info.memory_efficient_given ) {
			StrFloRow dm;
						//open infile

						// THE DATA WE WILL PROCESS
						std::vector<Sequence> seqs;
						std::vector<string> names;
						std::vector<DNA_b128_String> b128seqs;

						Extrainfos extrainfos;

						//for each dataset in the files


						for ( int ds = 0 ; ds < ndatasets || args_info.input_format_arg == input_format_arg_xml ; ds++ ){
							//no bootstrapping
							std::string runId("");
							if ( !no_incl_orig && numboot == 0){//only need to create one distance matrix
								if ( ! istream->read(b128seqs,runId,names,extrainfos)) {
									break;
								}

								const size_t numberOfSequences = b128seqs.size();
								ostream->printStartRun(names,runId,extrainfos);
								ostream->printHeader(numberOfSequences);

								for(size_t i = 0; i < numberOfSequences; ++i){
									fillMatrixRow(dm, b128seqs, trans_model, i, true);
									dm.setIdentifier(names.at(i));
									if(useFixFactor) applyFixFactorRow(dm,fixfactor);
									ostream->printRow(dm, names.at(i), i, true);
								}
							}
							//bootstrapping
							else{
								//TODO: implement bootstraping with the improved version of FastDist
								StrFloRow dm;
								//read original sequences
								if ( ! istream->readSequences(seqs,runId,extrainfos)) break;
								names.clear();names.reserve(seqs.size());
								for( size_t i=0;i<seqs.size();i++)
								{
									names.push_back(seqs[i].name);
								}
								const size_t numberOfSequences = seqs.size();
								ostream->printStartRun(names,runId,extrainfos);
								ostream->printHeader(numberOfSequences);
								if ( !no_incl_orig ){//create the distance matrix for the original sequences
									Sequences2DNA_b128(seqs,b128seqs);
									for(size_t i = 0; i < numberOfSequences; ++i){
										fillMatrixRow(dm, b128seqs, trans_model, i, true);
										dm.setIdentifier(names.at(i));
										if(useFixFactor) applyFixFactorRow(dm,fixfactor);
										ostream->printRow(dm, names.at(i), i, true);
									}
								}
								//start the bootstrapping
								//	  vector<Sequence> bootsequences;
								for ( int b = 0 ; b < numboot ; b++ ){
									bootstrapSequences(seqs,b128seqs);
									std::cout << numberOfSequences << std::endl;
									ostream->printBootstrapSpliter(numberOfSequences);
									for(size_t i = 0; i < numberOfSequences; ++i){
										fillMatrixRow(dm, b128seqs, trans_model, i, true);
										dm.setIdentifier(names.at(i));
										if(useFixFactor) applyFixFactorRow(dm,fixfactor);
										ostream->printRow(dm, names.at(i), i, true);
									}
								}
							}
							ostream->printEndRun();
						}//end data set loop
		}
		else{
			StrDblMatrix dm;
			//open infile
			// THE DATA WE WILL PROCESS
			std::vector<Sequence> seqs;
			std::vector<string> names;
			std::vector<DNA_b128_String> b128seqs;
			Extrainfos extrainfos;

			//for each dataset in the files
			for ( int ds = 0 ; ds < ndatasets || args_info.input_format_arg == input_format_arg_xml ; ds++ ){
				//no bootstrapping
				std::string runId("");
				if ( !no_incl_orig && numboot == 0){//only need to create one distance matrix
					if ( ! istream->read(b128seqs,runId,names,extrainfos)) break;
					fillMatrix(dm, b128seqs, trans_model);
					ostream->printStartRun(names,runId,extrainfos);
					//          freeXmlStrings(extrainfos);
					dm.setIdentifiers(names);
					if(useFixFactor) applyFixFactor(dm,fixfactor);
					//	  output << dm;
					ostream->print(dm);
				}
				//bootstrapping
				else{
					//read original sequences

					if ( ! istream->readSequences(seqs,runId,extrainfos)) break;

					names.clear();names.reserve(seqs.size());
					for( size_t i=0;i<seqs.size();i++)
						names.push_back(seqs[i].name);

					ostream->printStartRun(names,runId,extrainfos);
					//          freeXmlStrings(extrainfos);
					if ( !no_incl_orig ){//create the distance matrix for the original sequences
						Sequences2DNA_b128(seqs,b128seqs);
						fillMatrix(dm, b128seqs, trans_model);
						dm.setIdentifiers(names);
						if(useFixFactor)
							applyFixFactor(dm,fixfactor);
						//       	       output << dm;
						ostream->print(dm);
					}
					//start the bootstrapping
					//	  vector<Sequence> bootsequences;
					for ( int b = 0 ; b < numboot ; b++ ){
						//Sequence::bootstrapSequences(seqs,bootsequences);
						//Sequences2DNA_b128(bootsequences,b128seqs);
						bootstrapSequences(seqs,b128seqs);
						fillMatrix(dm, b128seqs, trans_model);
						dm.setIdentifiers(names);
						if(useFixFactor) applyFixFactor(dm,fixfactor);
						ostream->print(dm);
					}
				}
				ostream->printEndRun();
			}//end data set loop

			//OUTPUT THE TREES
		}////Mehmood's Changes End : email: malagori@kth.se
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

