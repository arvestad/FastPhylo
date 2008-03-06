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
#include "string_compare.hpp"

using namespace std;



void
print_options(char *note = NULL){
  if ( note != NULL ){
    cout << "ERROR: " << note <<endl<< endl;
  }
  cout <<
    "OPTIONS        MAN/DEF   DESCRIPTION\n"
    //    " -D val          K2P     Distance function 'val' should be either JC,K2P,FAKE_F84,TN93,HAMMING\n"
    " -i infiles       *      Phylip sequence files to produce matrices for\n"
    " -o file                 Out file for matrix/matrices. If not specified ends up in infile.dist\n"
"\n"
    " -h or --help            Print help message\n"
    "\n"
    " -b num                  Bootstrap num times and create matrix for each.\n"
    " --no-incl-orig          If the distance matrix from the original sequences should be included\n"
    " --seed num              Random seed. If not specified the current timestamp will be used.\n"
    "\n"
    " --tstvratio     2.0     What transition/transvertion ratio to use.\n"
    " --pyrtvratio    2.0     For the TN model there are two ratios. \'--tstvratio\' is the ratio for purine transitions while \'--pyrtvratio\' is for pyrimidines.\n"
    " --no-tstvratio          If given fixed ts/tv ratios will not be used.\n"
    "\n"
    " --fixfactor float       A float specifying what factor to use for saturated data. If not given -1 in the entry.\n"
    //" --combine               Combine all input files to one mutliple alignment set.\n"
    " -m ndatasets    1       Number of datasets in each file. (can not be used with --combine)\n"
  << endl;

  cout << "EXAMPLE USAGE" << endl;
  cout << " $./naivedist -i infile.seq" <<endl;
  cout << " $./naivedist -i infile.seq -o outfile.txt" <<endl;
  cout << " Distance matrices for all files matching *.seq are created" << endl; 
  cout << " $./naivedist -i *.seq -o outfile.txt" <<endl;
  cout << " Creates 10 bootstrapped matrices " <<endl; 
  cout << " $./naivedist -i *.seq -o outfile.txt -b 10 --no-incl-orig" <<endl;
  
  
  exit(1);
}

void 
fillMatrix_K2P_naive(StrDblMatrix &dm, std::vector<Sequence> &seqs,
	       sequence_translation_model trans_model);


int
main(int argc,
     char **argv){

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
      std::vector<std::string> names;

      //for each dataset in the files
      for ( int ds = 0 ; ds < ndatasets ; ds++ ){
	//read original sequences
	Sequence::readSequences(seqs,fin);
	names.clear();names.reserve(seqs.size());
	for( size_t i=0;i<seqs.size();i++)
	  names.push_back(seqs[i].name);
	
	if ( !no_incl_orig ){//create the distance matrix for the original sequences
	  fillMatrix_K2P_naive(dm, seqs, trans_model);
	   dm.setIdentifiers(names);
	   if(useFixFactor) applyFixFactor(dm,fixfactor);
	   printPHYLIPfast(dm,outfile);
	 }
	 //start the bootstrapping
	 vector<Sequence> bootsequences;
	 for ( int b = 0 ; b < numboot ; b++ ){
	   Sequence::bootstrapSequences(seqs,bootsequences);
	   fillMatrix_K2P_naive(dm, bootsequences, trans_model);
	  dm.setIdentifiers(names);
	  if(useFixFactor) applyFixFactor(dm,fixfactor);
	  printPHYLIPfast(dm,outfile);
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

  return 0;
}



//------------------------------------------------------------------
//FILL MATRIX


void 
fillMatrix_K2P_naive(StrDblMatrix &dm, std::vector<Sequence> &seqs,
	       sequence_translation_model trans_model){

  const size_t numSequences = seqs.size();
  const size_t strlen = seqs[0].seq.size();

  dm.resize(numSequences);
  
  for ( size_t i = 0 ; i < numSequences ; i++ ){
    dm.setDistance(i,i,0);
    ML_string_distance  ml_dist;
    for ( size_t j = i+1 ; j < numSequences ; j++ ){
      float divergence_matrix[4][4];
      int dels = complete_dna_string_compare(divergence_matrix,seqs[i].seq,seqs[j].seq);
      TN_string_distance tndist = divergence_matrix_2_TN_distance(divergence_matrix,dels);
      simple_string_distance sd = convert_TN_string_distance_to_simple(tndist);

      if ( trans_model.no_tstvratio )
        ml_dist = compute_K2P(strlen,sd);
      else
        ml_dist = compute_K2P_fixratio(strlen,sd,trans_model.tstvratio);
     
      dm.setDistance(i,j,ml_dist.distance);
    }
  }

}












