#include "fastprot_mpi_gengetopt.h"
#include "log_utils.hpp"
#include "file_utils.hpp"
#include "PhylipMaInputStream.hpp"
#include "DataInputStream.hpp"
#include "FastaInputStream.hpp"
#include "DataOutputStream.hpp"
#include "XmlOutputStream.hpp"
#include "fileFormatSchema.hpp"
#include "ProtDistCalc.hpp"
#include "ProtSeqUtils.hpp"
#include "../../DistanceMatrix.hpp"

#include <mpi.h>

#include <string>
#include <vector>
#include <cstdlib>
#include <numeric>

#include "mpi.h"

#ifdef WITH_LIBXML
#include "XmlInputStream.hpp"
#endif //WITH_LIBXML


int main (int argc, char **argv){
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
		gengetopt_args_info args_info;


#ifndef WITH_LIBXML
		if (args_info.input_format_arg == input_format_arg_xml){
			std::cerr << "The software was built with WITH_LIBXML=OFF. Please rebuild it if you want XML functionality." << std::endl; exit(EXIT_FAILURE);
		}
#endif // WITH_LIBXML

		if (cmdline_parser(argc, argv, &args_info) != 0)
			exit(EXIT_FAILURE);

		if (args_info.print_relaxng_input_given && args_info.print_relaxng_output_given){
			std::cerr << "error: --print-relaxng-input and --print-relaxng-output can not be used at the same time" << std::endl; exit(EXIT_FAILURE);
		}

		if (args_info.print_relaxng_input_given) {
			std::cout << fastphylo_prot_sequence_xml_relaxngstr << std::endl;
			exit(EXIT_SUCCESS);
		}
		if (args_info.print_relaxng_output_given) {
			std::cout << fastphylo_distance_matrix_xml_relaxngstr << std::endl;
			exit(EXIT_SUCCESS);
		}

		if (args_info.number_of_runs_given && args_info.input_format_arg != input_format_arg_phylip) {
			std::cerr << "error: --number-of-runs can only be used together with --input-format=phylip" << std::endl;
			exit(EXIT_FAILURE);
		}

		//--------------------------------------------------------------
		// Read translation model

		// prot_sequence_translation_model trans_model;

		if (! args_info.model_file_given){
			switch (args_info.distance_function_arg) {
			case distance_function_arg_ID : trans_model.model = id; break;
			case distance_function_arg_JC : trans_model.model = jc; break;
			case distance_function_arg_JCK : trans_model.model = jck; break;
			case distance_function_arg_JCSS : trans_model.model = jcss; break;
			case distance_function_arg_WAG : trans_model.model = wag; break;
			case distance_function_arg_JTT : trans_model.model = jtt; break;
			case distance_function_arg_DAY : trans_model.model = day; break;
			case distance_function_arg_ARVE : trans_model.model = arve; break;
			case distance_function_arg_MVR : trans_model.model = mvr; break;
			default: std::cerr << "error: model chosen not available" << std::endl; exit(EXIT_FAILURE);
			}
		} else {
			// read file from args_info.model_file_arg
		}



		if (args_info.maximum_likelihood_given &&
				(trans_model.model == id || trans_model.model == jc ||
						trans_model.model == jck || trans_model.model == jcss )) {
			std::cerr << "error: --maximum-likelihood can not be used with --distance-function=ID, JC, JCK or JCSS" << std::endl;
			exit(EXIT_FAILURE);
		}

		trans_model.step_size = args_info.speed_arg;

		if (args_info.pfam_given)
			trans_model.tp = norm;
		else
			trans_model.tp = flat;

		trans_model.ml = args_info.maximum_likelihood_given;
		bool remove_indels = args_info.remove_indels_given;
		int ndatasets = args_info.number_of_runs_arg;

		//----------------------------------------------
		// BOOTSTRAPPING
		int numboot = args_info.bootstraps_arg;
		bool no_incl_orig = args_info.no_incl_orig_given;

		if ( args_info.seed_given )
			srand((unsigned int )args_info.seed_arg);
		else
			srand((unsigned int)time(NULL));

		try {
			char * inputfilename = 0;
			char * outputfilename = 0;

			DataInputStream *istream;
			DataOutputStream *ostream;

			switch( args_info.inputs_num )
			{  case 0: break; /* inputfilename will be null and indicate stdin as input */
			case 1: inputfilename =  args_info.inputs[0]; break;
			default: std::cerr << "Error: you can at most specify one input filename" << std::endl;
			exit(EXIT_FAILURE);
			}

			if( args_info.outfile_given )
			{  outputfilename = args_info.outfile_arg;  }

			switch ( args_info.input_format_arg )
			{
			case input_format_arg_fasta: istream = new FastaInputStream(inputfilename);  break;
			case input_format_arg_phylip: istream = new PhylipMaInputStream(inputfilename);  break;
#ifdef WITH_LIBXML
			case input_format_arg_xml: istream = new XmlInputStream(inputfilename); break;
#endif // WITH_LIBXML
			default: exit(EXIT_FAILURE);
			}

			switch ( args_info.output_format_arg )
			{
			case output_format_arg_phylip: ostream = new PhylipDmOutputStream(outputfilename);  break;
			case output_format_arg_xml: ostream = new XmlOutputStream(outputfilename); break;
			default: exit(EXIT_FAILURE);
			}

			StrDblMatrix dm;

			std::vector<Sequence> seqs;
			std::vector<std::string> names;
			Extrainfos extrainfos;

			starttime = MPI::Wtime();


			//for each dataset in the file
			for ( int ds = 0 ; ds < ndatasets || args_info.input_format_arg == input_format_arg_xml ; ds++ ){
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

				ostream->printEndRun();
			}

			endtime = MPI::Wtime();

			delete ostream;
			delete istream;
			//printf("time: %lf\n", endtime-starttime);
		}
		catch(...){
			throw;
		}
		cmdline_parser_free(&args_info);
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
