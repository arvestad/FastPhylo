//--------------------------------------------------
//                                        
// File: drawtree.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: drawtree.cpp,v 1.5 2006/12/10 19:58:53 isaac Exp $                                 
//
//--------------------------------------------------




#include <string>
#include <time.h>
#include "arg_utils_ext.hpp"
#include "file_utils.hpp"
#include <iostream>

#include <fstream>
#include "SequenceTree.hpp"

using namespace std; 

void
print_options(char *note = NULL){
  if ( note != NULL ){
    cout << "ERROR: " << note <<endl<< endl;
  }
  cout <<
    "OPTIONS           MAN/DEF   DESCRIPTION\n"
    " -seqs file                  A sequence file.\n"
    " -d int             1       Number of datasets in file.\n"
    " -otree file                A file to print the phylip tree with internal node names set.\n"
    " -seqgen                    Set the names of internal nodes as seqgen would\n"
    " -h or --help               Print this help message\n"
       << endl;

  cout << "EXAMPLE USAGE" << endl;
  cout << " $./drawtree treefile.txt -seqs infile.seq" <<endl;
  
  
  exit(1);
}


int
main( int argc, char ** argv){

  TRY_EXCEPTION();
  
  if ( HAS_OPTION("-h") || HAS_OPTION("--help") )
    print_options();

  //--------
  if ( ! file_exists(argv[1]) )
    print_options("tree input file doesn't exist");
  ifstream inTreeFile;
  open_read_stream(argv[1], inTreeFile);
  
  //------
  
  bool seqgen = GET_BOOLEAN_OPTION_VAL("-seqgen");
  //--------------
  
  SequenceTree resultTree;
  SequenceTree::NodeVector nodes;
  //files
  vector<Sequence> seqs;
  int numdatasets = 1;
  char *dstr = GET_OPTION_VAL("-d");
  if ( dstr != NULL ) numdatasets = atoi(dstr);

  
  char *seqfile = GET_OPTION_VAL("-seqs");
  ifstream inSeqFile;
  if (seqfile != NULL )
    open_read_stream(seqfile,inSeqFile);

  char *otreefile = GET_OPTION_VAL("-otree");
  ofstream outTreeFile;
  if ( otreefile != NULL )
    open_write_stream(otreefile,outTreeFile);
  
  for ( int i = 0 ; i < numdatasets ; i++ ){
    cout << "**************"<< endl;
    cout << "Dataset " << (i+1) << endl;
    
    //read the tree
    resultTree = SequenceTree(inTreeFile);
    if ( seqgen ){
      nodes.clear();
      resultTree.addNodesInPrefixOrder(nodes);
      int nodename = resultTree.getNumLeafs() + 1;
      for ( size_t i = 0 ; i < nodes.size() ; i++ ){
        if ( nodes[i]->isLeaf() ) continue;        
        nodes[i]->data.s.name = string("")+nodename;
        nodename++;
      }
    }
    else{
      resultTree.setNodeNames();
    }
    if ( seqfile != NULL ) {
      resultTree.mapSequencesOntoTree(inSeqFile);
    }
    resultTree.drawTree(cout);
    if ( seqfile != NULL )
      resultTree.printSequences(cout);
    if ( otreefile != NULL )
      outTreeFile << resultTree << endl;

  }
  cout << "**************"<< endl;
  
  CATCH_EXCEPTION();

  return 1;
}









