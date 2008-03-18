///////////////////////////////////////////////
//                                           //
// File: fastdist.cpp                             //
//                                           //
// Created on: <14-Mar-2005 14:21:14 cherub>                            //
//                                           //
// Author: Isaac Elias                       //
// Email: isaac@nada.kth.se                  //
///////////////////////////////////////////////

#include "Sequences2DistanceMatrix.hpp"
#include <string>
#include <iostream>
#include <time.h>
#include <fstream>
#include <assert.h>

#include "file_utils.hpp"
#include <iomanip>
#include "log_utils.hpp"
#include "fastdist_gengetopt.h"
#include "NeighborJoining.hpp"
#include "DataInputStream.hpp"
#include "XmlInputStream.hpp"
#include "DataOutputStream.hpp"
#include "XmlOutputStream.hpp"

using namespace std;

int
main(int argc,
     char **argv){

  TRY_EXCEPTION();

  sequence_translation_model trans_model;
  gengetopt_args_info args_info;

  if (cmdline_parser (argc, argv, &args_info) != 0)
    exit(EXIT_FAILURE);

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
  int ndatasets = args_info.datasets_arg;

  //FINNISHED PARSING ARGS
  //---------------------------------------------------------
  // START BUILING MATRICES
  //

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
    case input_format_arg_phylip_multialignment: istream = new PhylipMaInputStream(inputfilename);  break;
    case input_format_arg_xml: istream = new XmlInputStream(inputfilename); break;
   default: exit(EXIT_FAILURE);
}

  switch ( args_info.output_format_arg )
    {
    case output_format_arg_phylip_dm: ostream = new PhylipDmOutputStream(outputfilename);  break;
    case output_format_arg_xml: ostream = new XmlOutputStream(outputfilename); break;
   default: exit(EXIT_FAILURE);
}


    StrDblMatrix dm;
      //open infile
  
       // THE DATA WE WILL PROCESS
      std::vector<Sequence> seqs;
      std::vector<string> names;
      std::vector<DNA_b128_String> b128seqs;

      //for each dataset in the files

      for ( int ds = 0 ; ds < ndatasets || args_info.input_format_arg == input_format_arg_xml ; ds++ ){
	//no bootstrapping
	if ( !no_incl_orig && numboot == 0){//only need to create one distance matrix
          if ( ! istream->read(names,b128seqs)) break;
	  fillMatrix(dm, b128seqs, trans_model);
	  dm.setIdentifiers(names);
	  if(useFixFactor) applyFixFactor(dm,fixfactor);
	  //	  output << dm;
	  ostream->print(dm);

	}
	//bootstrapping
	else{
	  //read original sequences

          if ( ! istream->readSequences(seqs)) break;

	  names.clear();names.reserve(seqs.size());
	  for( size_t i=0;i<seqs.size();i++)
	    names.push_back(seqs[i].name);

  
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
      }//end data set loop

  //OUTPUT THE TREES


      delete ostream;
      delete istream;

  }



  catch(...){


    throw;
  }

  CATCH_EXCEPTION();
  return 0;
}

