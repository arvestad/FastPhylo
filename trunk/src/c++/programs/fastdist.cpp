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
#include "arg_utils_ext.hpp"
#include "file_utils.hpp"
#include <iomanip>
#include "log_utils.hpp"

using namespace std;



void
print_options(char *note = NULL){
  if ( note != NULL ){
    cout << "ERROR: " << note <<endl<< endl;
  }
  cout <<
    "OPTIONS        MAN/DEF   DESCRIPTION\n"
    " -D val          K2P     Distance function 'val' should be either JC,K2P,FAKE_F84,TN93,HAMMING\n"
    " -i infiles       *      Phylip sequence files to produce matrices for\n"
    " -o file                 Out file for matrix/matrices. If not specified ends up in infile.dist\n"
"\n"
    " -h or --help            Print help message\n"
    "\n"
    " -b num                  Bootstrap num times and create matrix for each.\n"
    " --no-incl-orig          If the distance matrix from the original sequences should be included\n"
    " --seed num              Random seed. If not specified the current timestamp will be used.\n"
    "\n"
    " --no-ambiguities        Ignore ambiguities.\n"
    " --no-ambig-resolve      Specifies that ambigious symbols should not be resolved by nearest neighbor.\n"
    " --no-transprob          Specifies that the transition probabilities should not be used in the ambiguity model\n"
    " --af val        UNI     Ambiguity frequency model. Either UNI for uniform or BASE for basefrequencies.\n"
    "\n"
    " --tstvratio     2.0     What transition/transvertion ratio to use.\n"
    " --pyrtvratio    2.0     For the TN model there are two ratios. \'--tstvratio\' is the ratio for purine transitions while \'--pyrtvratio\' is for pyrimidines.\n"
    " --no-tstvratio          If given fixed ts/tv ratios will not be used.\n"
    "\n"
    " --fixfactor float       A float specifying what factor to use for saturated data. If not given -1 in the entry.\n"
    //" --combine               Combine all input files to one mutliple alignment set.\n"
    " -m ndatasets    1       Number of datasets in each file. (can not be used with --combine)\n"
    "\n----------------------\nTODO\n"
    " --bootmode mode ?????     Kind of bootstrapping.\n";
  cout << endl;

  cout << "EXAMPLE USAGE" << endl;
  cout << " $./fastdist -i infile.seq" <<endl;
  cout << " $./fastdist -D JC -i infile.seq -o outfile.txt" <<endl;
  cout << " Distance matrices for all files matching *.seq are created" << endl; 
  cout << " $./fastdist -D JC -i *.seq -o outfile.txt" <<endl;
  cout << " Creates 10 bootstrapped matrices " <<endl; 
  cout << " $./fastdist -i *.seq -o outfile.txt -b 10 --no-incl-orig" <<endl;
  
  
  exit(1);
}

int
main(int argc,
     char **argv){

  TRY_EXCEPTION();
  //-----
  if ( HAS_OPTION("-h") || HAS_OPTION("--help") )
    print_options();


  //---

  sequence_translation_model trans_model;

  //-----------------------------------------------------
  // EVOLUTIONARY MODEL
  trans_model.model = K2P;
  char *df_str = GET_OPTION_VAL("-D");
  if ( df_str != NULL ){
    if ( STREQ("JC", df_str) )
      trans_model.model  = JC;
    else if ( STREQ("K2P",df_str) )
      trans_model.model = K2P;
    else if ( STREQ("TN93",df_str) )
      trans_model.model = TN93;
    else if ( STREQ("HAMMING",df_str) )
      trans_model.model = HAMMING_DISTANCE;
    else
      print_options("Unkown option to -D");
  }
  
  trans_model.no_tstvratio = false;
  if ( HAS_OPTION("--no-tstvratio") ){
    trans_model.no_tstvratio = true;
  }

  trans_model.tstvratio = 2.0;
  SET_FLOAT_OPTION_VAL("--tstvratio",&(trans_model.tstvratio));
 
  trans_model.pyrtvratio = 2.0;
  SET_FLOAT_OPTION_VAL("--pyrtvratio",&(trans_model.pyrtvratio));
  
  //----------------------------------------------
  // BOOTSTRAPPING
  int numboot = 0;
  SET_INT_OPTION_VAL("-b",&numboot);

  bool no_incl_orig = false;
  if ( HAS_OPTION("--no-incl-orig") ){
    no_incl_orig = true;
  }

  char *seedval = GET_OPTION_VAL("--seed");
  if ( seedval == NULL )
    srand((unsigned int)time(NULL));
  else
    srand((unsigned int)atoi(seedval));
  //-----------------------------------------------
  // AMBIGUITIES
  trans_model.no_ambiguities = false;
  if ( HAS_OPTION("--no-ambiguities") ){
    trans_model.no_ambiguities = true;
  }
  trans_model.no_ambig_resolve = false;
  if ( HAS_OPTION("--no-ambig-resolve") ){
    trans_model.no_ambig_resolve = true;
  }

  
  trans_model.no_transition_probs = false;
  if ( HAS_OPTION("--no-transprob") ){
    trans_model.no_transition_probs = true;
  }

  char *af_str = GET_OPTION_VAL("--af");
  trans_model.use_base_freqs = false;
  if ( af_str != NULL ){
    if ( STREQ("UNI", af_str) )
      trans_model.use_base_freqs = false;
    else if ( STREQ("BASE",af_str) )
      trans_model.use_base_freqs = true;
    else
      print_options("Unkown option to --af");
  }
  
  //----------------------------------------------
  bool useFixFactor = false;
  float fixfactor=1;
  char *ffstr = GET_OPTION_VAL("--fixfactor");
  if ( ffstr != NULL ){
    fixfactor = atof(ffstr);
    useFixFactor = true;
  }
  //---

  int ndatasets = 1;
  SET_INT_OPTION_VAL("-m",&ndatasets);
 
  //----------------------------------------------
  // SEQUENCE FILES
  vector<char *> infiles;
  if ( !GET_LIST_OF_ARGS("-i",infiles) )
    print_options("No input file");

  char *outfile_name = GET_OPTION_VAL("-o");
 
  //FINNISHED PARSING ARGS
  //---------------------------------------------------------
  // START BUILING MATRICES
  //
  ifstream fin;
  FILE *outfile = NULL;
  if( outfile_name!=NULL )//only one output file for all 
    outfile = open_write_file(outfile_name);
  
  try {
    StrDblMatrix dm;
  
    for ( size_t i=0; i<infiles.size(); i++ ){
      //open infile
      open_read_stream(infiles[i],fin);
      //if no name for output file was used open infile+".dist"
      if ( outfile_name==NULL ){
	string outname(infiles[i]);
	outname += ".dist";
	outfile = open_write_file(outname.c_str());
      }
    
      // THE DATA WE WILL PROCESS
      std::vector<Sequence> seqs;
      std::vector<string> names;
      std::vector<DNA_b128_String> b128seqs;

      //for each dataset in the files
      for ( int ds = 0 ; ds < ndatasets ; ds++ ){
	//no bootstrapping
	if ( !no_incl_orig && numboot == 0){//only need to create one distance matrix
	  DNA_b128_StringsFromPHYLIP(fin,names,b128seqs);
	  fillMatrix(dm, b128seqs, trans_model);
	  dm.setIdentifiers(names);
	  if(useFixFactor) applyFixFactor(dm,fixfactor);
	  printPHYLIPfast(dm,outfile);
	}
	//bootstrapping
	else{
	  //read original sequences
	  Sequence::readSequences(seqs,fin);
	  names.clear();names.reserve(seqs.size());
	  for( size_t i=0;i<seqs.size();i++)
	    names.push_back(seqs[i].name);
	  
	  if ( !no_incl_orig ){//create the distance matrix for the original sequences
	    Sequences2DNA_b128(seqs,b128seqs);
	    fillMatrix(dm, b128seqs, trans_model);
	    dm.setIdentifiers(names);
	    if(useFixFactor) applyFixFactor(dm,fixfactor);
	    printPHYLIPfast(dm,outfile);
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
	    printPHYLIPfast(dm,outfile);
	  }
	}
      }//end data set loop
      
      // CLOSE FILES so that they can be reused in the next interation
      fin.close();
      if( outfile_name==NULL ) fclose(outfile);
    }//end infile loop
  }
  catch(...){
    fclose(outfile);
    throw;
  }

  if ( outfile_name!=NULL )
    fclose(outfile);

  CATCH_EXCEPTION();
  return 0;
}












