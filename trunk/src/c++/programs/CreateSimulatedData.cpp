

//--------------------------------------------------
//                                        
// File: Clustal2gaplessPhylip.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: CreateSimulatedData.cpp,v 1.10 2006/12/05 09:08:29 isaac Exp $                                 
//
//--------------------------------------------------



#include <string>
#include <time.h>
#include "arg_utils_ext.hpp"
#include "file_utils.hpp"
#include <iostream>
#include <vector>
#include "Simulator.hpp"
#include "SequenceTree.hpp"

#include <fstream>
#include "stl_utils.hpp"

using namespace std; 

void
print_options(char *note = NULL){
  if ( note != NULL ){
    cout << "ERROR: " << note <<endl<< endl;
  }
  cout <<
    "OPTIONS                MAN/DEF   DESCRIPTION\n"
    " --leafs int int int       10     The different leaf sizes (list of ints)\n"
    " --diamf flt flt flt     0.25    Diameter factor (list of flts)\n"
    " --seqlen int int int     200    Sequence length\n"
    " --numtrees int            20    The number of trees for each size\n"
    " --wa                            Flag to write ancestral sequences\n"
    " --contract                      Contracts all internal edges for which the number of events was 0\n"
    " --prefix str              \"\"    The file prefix for all output files\n"  
    " --unrooted                      Flag that the outputed tree should be unrooted\n"
  
    "\n -h or --help               Print this help message\n"
    "\n--------------------\nTODO:\n --ultrametricDeviation\n";
  cout << endl;

  cout << "EXAMPLE USAGE" << endl;
  cout << " $./CreateSimulatedData --leafs 10 20 30 --diamf 0.05 0.2 --seqlen 100 200 --numtrees 2" <<endl;
  
  
  exit(1);
}




int
main( int argc, char ** argv){

  if  ( HAS_OPTION("-h") || HAS_OPTION("--help") )
    print_options();

  
  vector<int> leafSizes;
  if( !HAS_OPTION("--leafs") ){
    leafSizes.push_back(10);
  }
  else{
    GET_LIST_OF_INTS("--leafs", leafSizes);
  }

  vector<float> diameterFactors;
  if( !HAS_OPTION("--diamf") ){
    diameterFactors.push_back(0.25);
  }
  else{
    GET_LIST_OF_FLOATS("--diamf", diameterFactors);
  }


  vector<int> seqLens;
  if( !HAS_OPTION("--seqlen") ){
    seqLens.push_back(200);
  }
  else{
    GET_LIST_OF_INTS("--seqlen", seqLens);
  }  

  int numTrees = 20;
  char *tmp = GET_OPTION_VAL("--numtrees");
  if( tmp!=NULL )
    numTrees = atoi(tmp);
  
  bool writeAncestors = GET_BOOLEAN_OPTION_VAL("--wa");
  bool contract = GET_BOOLEAN_OPTION_VAL("--contract");
  bool unrooted = GET_BOOLEAN_OPTION_VAL("--unrooted");

  string prefix;
  tmp = GET_OPTION_VAL("--prefix");
  if( tmp!=NULL )
    prefix.append(tmp);

  int totalNumberOfCases = leafSizes.size()*seqLens.size()*diameterFactors.size()*numTrees;
  int currentCase =0;

  for(size_t leafi=0 ; leafi<leafSizes.size(); leafi++){
    for( size_t seqi=0 ; seqi<seqLens.size() ; seqi++){
      for( size_t diami=0 ; diami<diameterFactors.size() ; diami++){
	string treeFileName = prefix;
	treeFileName +="tree_";
	createFileName(treeFileName, leafSizes[leafi],4, seqLens[seqi], diameterFactors[diami]);
	ofstream treeFile;
	open_write_stream(treeFileName.c_str(),treeFile);
	
	string seqFileName = prefix;
	seqFileName +="sequences_";
	createFileName(seqFileName, leafSizes[leafi],4, seqLens[seqi], diameterFactors[diami]);
	ofstream seqFile;
	open_write_stream(seqFileName.c_str(),seqFile);
	
	for( int treenum=0 ; treenum< numTrees ; treenum++){
	  SEPARATOR();
	  PRINT(currentCase++);PRINT(totalNumberOfCases);
	  PRINT(leafSizes[leafi]);PRINT(seqLens[seqi]);
	  PRINT(diameterFactors[diami]);PRINT(treenum);
	  PRINT(treeFileName);

	  SequenceTree tree;
	  PRINT(tree);
	  LINE();
	  createRandomTreeUsingBeep(tree,leafSizes[leafi], 4);
	  LINE();
	  PRINT(tree);
	  LINE();
	  evolveSequencesUsingSeqGen(tree, seqLens[seqi], diameterFactors[diami],true);
	  if(contract){
	    tree.computeEdgeLengths();
	    tree.drawTree(cout);
	    int numContracted = tree.contractEdgesShorterThan(0);
	    PRINT(numContracted);
	    tree.drawTree(cout);
	  }
	  if(unrooted){
	    tree.shortcutDegree2Nodes();
	  }
	  SequenceTree::NodeVector nodes;
	  if(writeAncestors)
	    tree.addNodesInPrefixOrder(nodes);
	  else{
	    tree.addNodesInPrefixOrder(nodes);
	    for(size_t nodei=0;nodei<nodes.size();nodei++)
	      if(!nodes[nodei]->isLeaf())
		NAME(nodes[nodei])="";
	    nodes.clear();
	    tree.addLeafs(nodes);
	  }
	  LINE();
	  treeFile << tree<< endl;
	  LINE();
	  SequenceTree::printSequencesPhylip(nodes,seqFile);
	  LINE();
	}
	treeFile.close();
	seqFile.close();
      }
    }
    
  }
  
  
  return 1;
}
