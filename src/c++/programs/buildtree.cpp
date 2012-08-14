//--------------------------------------------------
//                                        
// File: buildtree.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: buildtree.cpp,v 1.19 2006/12/26 13:42:28 isaac Exp $                                 
//
//--------------------------------------------------

#include <string>
#include <time.h>
#include "arg_utils_ext.hpp"
#include "file_utils.hpp"
#include <iostream>

#include "Sequences2DistanceMatrix.hpp"

#include <string>
#include <iostream>
#include <time.h>
#include <fstream>
#include "NeighborJoining.hpp"
#include "stl_utils.hpp"
#include "log_utils.hpp"

typedef __gnu_cxx::hash_map<const SequenceTree , int, objhash, objeq> tree2int_map;

using namespace std; 



//
// Builds trees using the supplied distance methods.
// Each created tree is added to the tree2count map.
//
void
buildTrees(StrDblMatrix &dm, tree2int_map &tree2count, std::vector<NJ_method> &methods, str2int_hashmap &name2id){

  SequenceTree tree;
  for(size_t i=0 ; i<methods.size() ; i++){
    computeNJTree(dm,tree,methods[i]);
    tree.makeCanonical(name2id);
    tree2int_map::iterator iter = tree2count.find(tree);
    if(iter!=tree2count.end())
      (*iter).second++;
    else
      tree2count[tree] = 1;
  }


}


void
print_options(char *note = NULL){
  if ( note != NULL ){
    cout << "ERROR: " << note <<endl<< endl;
  }
  cout <<
    "OPTIONS           MAN/DEF     DESCRIPTION\n"
    "-mats <file>      or -seqs    Input file with distance matrices.\n"
    "-seqs <file>      or -mats    Sequence file.\n"
    //    "-d <integer>        1         Number of datasets in the file.\n"
    "-o <file>         %.tree      Output files. \n"
    "-no-counts                    If tree file should not contain the counts\n"
    "                              for the tree (PHYLIP compliance).\n"
    "\n"
    "-methods <strings>   NJ       The reconstruction methods to apply {NJ,FNJ,BIONJ}\n"
    "\n\nOptions for sequence files only.\n"
    "-boot <int>          0        Number of times to bootstrap.\n"
    "-seed <int>         time      Seed for random numbers.\n"
    "-exclude-orig                 Flag to exclude the original sequences set.\n"
    "-model <string>     K2P       Sequence evolution model to use for distance estimation {HAMMING,JC,K2P,TN93.\n"
    " -tstv              2.0       The transition/transvertion ratio to use.\n"
    " -pyrtv             2.0       For the TN model there are two ratios. \'-tstv\' is the ratio for purine transitions while \'--pyrtv\' is for pyrimidines.\n"
    "-fix-factor <int>    2        Distance for over saturated data is set to <int>*maxRealDistance.\n"
    "\n"
    "-h or --help Print this help message.\n"
    "\n\n"
    "EXAMPLE USAGE\n" 
    "The matirx in the input file is used to build three trees using the respective methods.\n"
    "$./buildtree -mats infile -d 20 -methods NJ BIONJ FNJ\n"
    "\nCompute distance matrix using TN93,apply the methods, do an additional 100 bootstraps\n"
    "and only output the unique topologies.\n"
    "$./buildtree -seqs infile -methods NJ BIONJ FNJ -boot 100 -model TN93 -unique-trees\n"
       << endl;  

  exit(1);
}


int
main(int argc, char **argv){
    
  if ( HAS_OPTION("-h") || HAS_OPTION("--help") )
    print_options();
  
  //------------------------------------------
  // FILES
  bool isMatrixFile = true;
  
  char *infileName = GET_OPTION_VAL("-mats");
  if( infileName==NULL ){
    infileName = GET_OPTION_VAL("-seqs");
    isMatrixFile = false;
  }
  if( infileName==NULL )
    print_options("No input file with matrices or sequences.\n");
  if( !file_exists(infileName) )
    print_options("Infile doesn't exist.");

  std::string outfileName(infileName);
  outfileName += ".tree";
  char *outarg = GET_OPTION_VAL("-o");
  if(outarg!=NULL)
    outfileName = std::string(outarg);

  //----------------------------------------------
  // DISTANCE METHODS
  std::vector<NJ_method> methods;
  std::vector<char*> methodNames;
  if(GET_LIST_OF_ARGS("-methods",methodNames)==false)
    methods.push_back(NJ);//defaul

  for( size_t i=0 ; i<methodNames.size() ; i++)
    if( STREQ("NJ", methodNames[i]) )
      methods.push_back(NJ);
    else if( STREQ("FNJ",methodNames[i]) )
      methods.push_back(FNJ);
    else if( STREQ("BIONJ",methodNames[i]) )
      methods.push_back(BIONJ);
    else
      print_options("Unkown option to -methods");

  //--------------
  bool noCounts = GET_BOOLEAN_OPTION_VAL("-no-counts");  
  //-----------------------------------------------
  // SEQUENCE MODEL
  sequence_translation_model trans_model;
   
  char *modelName = GET_OPTION_VAL("-model");
  if( modelName==NULL )
    modelName = "K2P";
 
  if ( STREQ("JC", modelName) )
    trans_model.model  = JC;
  else if ( STREQ("K2P",modelName) )
    trans_model.model = K2P;
  else if ( STREQ("TN93",modelName) )
    trans_model.model = TN93;
  else if ( STREQ("HAMMING",modelName) )
    trans_model.model = HAMMING_DISTANCE;
  else
    print_options("Unkown option to -model");

  //ts tv ratios
  trans_model.no_tstvratio = false; 
  trans_model.tstvratio = 2.0;
  SET_FLOAT_OPTION_VAL("-tstv",&(trans_model.tstvratio));
  trans_model.pyrtvratio = 2.0;
  SET_FLOAT_OPTION_VAL("--pyrtv",&(trans_model.pyrtvratio));
      
  float fixfactor=2;
  SET_FLOAT_OPTION_VAL("-fix-factor",&fixfactor);
  trans_model.no_ambiguities = false;
  if ( HAS_OPTION("--no-ambiguities") ){
    trans_model.no_ambiguities = true;
  }

  //ambiguity model
  trans_model.no_ambig_resolve = false;
  trans_model.no_transition_probs = false;
  trans_model.use_base_freqs = false;

  //------------------------------------
  // BOOTSTRAP FOR SEQUENCES
  int numBoot = 0;
  SET_INT_OPTION_VAL("-boot",&numBoot);
  
  bool excludeOriginal = GET_BOOLEAN_OPTION_VAL("-exclude-orig");
  
  char *seedval = GET_OPTION_VAL("-seed");
  if ( seedval == NULL )
    srand((unsigned int)time(NULL));
  else
    srand((unsigned int)atoi(seedval));
  
  //All options read
  //-----------------------------------------
  ifstream fin;
  ofstream fout;
  open_read_stream(infileName,fin);
  
  //create tree set.
  tree2int_map tree2count((size_t)(numBoot*1.3));

  str2int_hashmap name2id;
  std::vector<std::string> names;
  std::vector<DNA_b128_String> b128seqs;
  std::vector<Sequence> seqs;
  StrDblMatrix dm;

  if(!isMatrixFile){
    //if only output tree for actual sequences.
    if ( !excludeOriginal && numBoot==0){	  
      DNA_b128_StringsFromPHYLIP(fin,names,b128seqs);
      for(size_t namei=0 ; namei<names.size() ; namei++ )
	name2id[names[namei]] = namei;

      fillMatrix(dm, b128seqs, trans_model);
      dm.setIdentifiers(names);
      applyFixFactor(dm,fixfactor);
	  
      buildTrees(dm, tree2count, methods,name2id);
    }
    //bootstrap
    else{
      //read sequences
      Sequence::readSequences(seqs,fin);
      for(size_t namei=0 ; namei<seqs.size() ; namei++ ){
	name2id[seqs[namei].name] = namei;
	names.push_back(seqs[namei].name);
      }
	  
      if( !excludeOriginal ){
	Sequences2DNA_b128(seqs,b128seqs);
	fillMatrix(dm, b128seqs, trans_model);
	dm.setIdentifiers(names);
	applyFixFactor(dm,fixfactor);
	    
	buildTrees(dm, tree2count, methods,name2id);
      }
	 
      for(int b=0 ; b<numBoot ; b++ ){
	//boot and create dm.
	//	Sequence::bootstrapSequences(seqs,bootseqs);
	//Sequences2DNA_b128(bootseqs,b128seqs);
	bootstrapSequences(seqs,b128seqs);
	    
	fillMatrix(dm, b128seqs, trans_model);
	dm.setIdentifiers(names);
	applyFixFactor(dm,fixfactor);
	    
	buildTrees(dm, tree2count, methods,name2id);
	    
      }
    }
	
  }
  //IF the input file contains distance matrices
  else{
    dm.objInitFromStream(fin);
    for(size_t namei=0 ; namei<dm.getSize() ; namei++ )
      name2id[dm.getIdentifier(namei)] = namei;
    buildTrees(dm, tree2count, methods,name2id);

	
  }
  //OUTPUT THE TREES
  open_write_stream(outfileName.c_str(),fout);
  tree2int_map::iterator iter = tree2count.begin();
  for( ; iter!=tree2count.end() ; ++iter){
    if(noCounts)
      fout << (*iter).first << endl;
    else
      fout << (*iter).second << "  " << (*iter).first << endl; 
  }



      
  return 0;
}












