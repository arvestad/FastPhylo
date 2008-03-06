//--------------------------------------------------
//                                        
// File: sequence_nj.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: sequence_nj.cpp,v 1.4 2006/12/10 19:58:53 isaac Exp $                                 
//
//--------------------------------------------------

#include <string>
#include <time.h>
#include "arg_utils_ext.hpp"
#include "file_utils.hpp"
#include <iostream>

#include "SequenceBasedNJ.hpp"
#include <fstream>

using namespace std; 

void
print_options(char *note = NULL){
  if ( note != NULL ){
    cout << "ERROR: " << note <<endl<< endl;
  }
  cout <<
    "OPTIONS           MAN/DEF   DESCRIPTION\n"
    " -seqs file         *      A sequence file for which a tree should be built.\n"
    " -otree file               A name of the output tree file. Default stdout\n"
    " -d int             1      Number of datasets\n"               
    "\n"
    " -h or --help               Print this help message\n"
    " -V                         Print version\n";
  cout << endl;

  cout << "EXAMPLE USAGE" << endl;
  cout << " $./sequence_nj -seqs infile.seq" <<endl;
  
  
  exit(1);
}


int
main( int argc, char ** argv){

  TRY_EXCEPTION();
  
  if ( HAS_OPTION("-h") || HAS_OPTION("--help") )
    print_options();

  if ( HAS_OPTION("-V") )
    PRINTVERSION();


  //------
  char *seqfile = GET_OPTION_VAL("-seqs");
  if ( seqfile == NULL )print_options();

  
  int numdatasets = 1;
  char *dstr = GET_OPTION_VAL("-d");
  if ( dstr != NULL ) numdatasets = atoi(dstr);

  //seqs
  ifstream inSeqFile;
  open_read_stream(seqfile,inSeqFile);
  vector<Sequence> seqs;

  char *outfile_name = GET_OPTION_VAL("-otree");
  ofstream outTreeFile;
  if ( outfile_name != NULL )
    open_write_stream(outfile_name,outTreeFile);
  

  for ( int i = 0 ; i < numdatasets ; i++ ){
    seqs.clear();
    Sequence::readSequences(seqs,inSeqFile);
    SequenceTree result;
    computeSequenceBasedNJ(seqs,result);
    if ( outfile_name != NULL )
      outTreeFile << result << endl;
    else
      cout << result << endl;
  }
  //------------------

  if ( outfile_name != NULL )
    outTreeFile.close();
  inSeqFile.close();

  CATCH_EXCEPTION();

  return 1;
}






