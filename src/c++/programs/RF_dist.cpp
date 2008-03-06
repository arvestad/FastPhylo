//--------------------------------------------------
//                                        
// File: drawtree.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: RF_dist.cpp,v 1.3 2006/12/10 19:58:53 isaac Exp $                                 
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
    " -h or --help               Print this help message\n"
    " -d int              1      The number of trees in each file\n";
  cout << endl;

  cout << "EXAMPLE USAGE" << endl;
  cout << " $./RF_dist tree1.txt tree2.txt" <<endl;
  cout << " $./RF_dist tree1.txt tree2.txt -d 20" <<endl;
  
  
  exit(1);
}


int
main( int argc, char ** argv){

  TRY_EXCEPTION();
  
  if ( HAS_OPTION("-h") || HAS_OPTION("--help"))
    print_options();
  if ( argc < 3 )
    print_options("To few arguments");

  //--------
  if ( ! file_exists(argv[1]) )
    print_options("tree input file doesn't exist arg=1");
  ifstream inTreeFile1;
  open_read_stream(argv[1], inTreeFile1);

  if ( ! file_exists(argv[2]) )
    print_options("tree input file doesn't exist arg=2");
  ifstream inTreeFile2;
  open_read_stream(argv[2], inTreeFile2);
  
  
  SequenceTree tree1;
  SequenceTree tree2;

  SequenceTree::NodeVector nodes;
  
  int numdatasets = 1;
  char *dstr = GET_OPTION_VAL("-d");
  if ( dstr != NULL ) numdatasets = atoi(dstr);

    
  for ( int i = 0 ; i < numdatasets ; i++ ){
    cout << "**************"<< endl;
    cout << "Dataset " << (i+1) << endl;
    
    //read the trees
    tree1 = SequenceTree(inTreeFile1);
    cout << tree1 << endl;
    tree2 = SequenceTree(inTreeFile2);
    cout << tree2 << endl;
    cout << "RF dist " << SequenceTree::computeRobinsonFoulds(tree1,tree2) << endl;
  }
  cout << "**************"<< endl;
  

  CATCH_EXCEPTION();

  return 1;
}














