//--------------------------------------------------
//                                        
// File: aml_tree.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: aml_tree.cpp,v 1.9 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------


#include <string>
#include <time.h>
#include "arg_utils_ext.hpp"
#include "file_utils.hpp"
#include <iostream>

#include "dna_pairwise_sequence_likelihood.hpp"
#include "SequenceTree.hpp"
#include "Big_AML.hpp"
#include <fstream>
#include "AML_LeafLifting.hpp"
#include "AML_given_edge_probabilities.hpp"
#include "AML_local_improve.hpp"

using namespace std; 

void
print_options(char *note = NULL){
  if ( note != NULL ){
    cout << "ERROR: " << note <<endl<< endl;
  }
  cout <<
    "OPTIONS           MAN/DEF   DESCRIPTION\n"
    " -tree file                 A file containing a phylogenetic tree in phylip format. If not specified one is built.\n"
    " -seqs files         *      One or more sequence files which maps names in the tree onto sequences or for which a tree should be built.\n"
    " -otree file                A name of the output tree file. Default the input seq file +\".amltree\"\n"
    " -oseqs file                A name of the output sequence file for the aml sequenences. Default the input seq file +\".amlseqs\"\n"
    " -model m           JC      Which sequence model to use P_DIST, JC\n"
    " -use-edges                 Computes optimal AML sequences by using the given edge lengths\n"
    " -seqgen                    Set internal node names as Seq-Gen would\n"
    " -no-local-improve          Don't perform local improvement heuristic.\n"
    " -use-ancestral            The input tree already has ancestral sequences so just compute the likelihood and do improvement\n"
    " -use-parsimony             Compute the parsimony solution and do local improvement\n"
    " -d int             1       Number of datasets in file.\n"
    " -v                         If should give verbose output.\n"
    "\n"
    " -h or --help               Print this help message\n"
       << endl;

  cout << "EXAMPLE USAGE" << endl;
  cout << " $./aml_tree -seqs infile.seq" <<endl;
  cout << " $./aml_tree -tree species.tree -seqs infile.seq" <<endl;
  
  
  exit(1);
}


int
main( int argc, char ** argv){

  TRY_EXCEPTION();
  
  if ( HAS_OPTION("-h") || HAS_OPTION("--help") )
    print_options();

  //------
  if ( ! HAS_OPTION("-seqs") )
    print_options();
  
  vector<char *> seqfiles;
  GET_LIST_OF_ARGS("-seqs", seqfiles);

  //-----
  sequence_model model = JC;
  char *df_str = GET_OPTION_VAL("-model");
  if ( df_str != NULL ){
    if ( STREQ("JC", df_str) )
      model = JC;
    else if ( STREQ("P_DIST",df_str) )
      model = P_DISTANCE;
    else
      print_options("Unkown option to -model");
  }

  //--------------
  
  ofstream outTreeFile;
  string outname(seqfiles[0]);
  outname += ".amltree";
  char *outfile_name = GET_OPTION_VAL("-otree");
  if ( outfile_name != NULL )
    outname = string(outfile_name);
  
  open_write_stream(outname.c_str(),outTreeFile);
  //-----------------
  
  ofstream outSeqFile;
  outname.clear();
  outname = string(seqfiles[0]);
  outname += ".amlseqs";
  outfile_name = GET_OPTION_VAL("-oseqs");
  if ( outfile_name != NULL )
    outname = string(outfile_name);
  
  open_write_stream(outname.c_str(),outSeqFile);
  //----------------

  bool use_edges = GET_BOOLEAN_OPTION_VAL("-use-edges");
  bool no_local_improve = GET_BOOLEAN_OPTION_VAL("-no-local-improve");
  bool seqgen = GET_BOOLEAN_OPTION_VAL("-seqgen");
  bool use_parsimony = GET_BOOLEAN_OPTION_VAL("-use-parsimony");
  bool use_ancestral = GET_BOOLEAN_OPTION_VAL("-use-ancestral");
  bool verbose = GET_BOOLEAN_OPTION_VAL("-v");
  //Options read
  //-------------------------------------------
  //start the algorithm
  SequenceTree resultTree;

  //files
  bool small_aml = HAS_OPTION("-tree");
  vector<Sequence> seqs;
  int numdatasets = 1;
  char *dstr = GET_OPTION_VAL("-d");
  if ( dstr != NULL ) numdatasets = atoi(dstr);
  
  ifstream inSeqFile;
  open_read_stream(seqfiles[0],inSeqFile);
  if ( seqfiles.size() > 1 ) PROG_ERROR("Not implmented seqfiles>1");
  ifstream inTreeFile;
  if (small_aml)
    open_read_stream(GET_OPTION_VAL("-tree"),inTreeFile);


  for ( int i = 0 ; i < numdatasets ; i++ ){
    cout << "**************"<< endl;
    cout << "Dataset " << (i+1) << endl;
    
    //If a tree is given as input 
    if ( small_aml ){
      //read the tree
      resultTree = SequenceTree(inTreeFile);
      //seqgen names for internal nodes
      if ( seqgen ){
        SequenceTree::NodeVector nodes;
        resultTree.addNodesInPrefixOrderLeftRight(nodes);
        int nodename = resultTree.getNumLeafs() + 1;
        for ( size_t i = 0 ; i < nodes.size() ; i++ ){
          if ( nodes[i]->isLeaf() ) continue;        
          nodes[i]->data.s.name = string("")+nodename;
          nodename++;
        }
      }
      
      cout << "shortcutting degree 2 nodes" << endl;
      resultTree.shortcutDegree2Nodes();
      resultTree.mapSequencesOntoTree(inSeqFile);
      
      if ( use_ancestral ){
        cout << "using ancestral sequences" <<endl;
      }
      else if ( use_parsimony ){
        cout << "using parsimony sequences" <<endl;
        resultTree.computeMostParsimoniousSequences();
      }
      else if ( use_edges ){
        cout << "using edge lengths in tree" <<endl;
        convert_edge_lengths_to_probabilities(resultTree,model);
        computeAML_given_edge_probabilities(resultTree);
      }
      else {
        cout << "using leaflifting" << endl;
        computeOptimal_AML_LeafLifting(resultTree, model);
      }
      double logl = resultTree.compute_loglikelihood();
      if ( verbose )
        resultTree.verbosePrint(cout);
      if ( ! no_local_improve ){
        cout <<"Before improvement loglikelihood = " << logl << endl;
        cout << "Doing local improvement heuristic" << endl;
        AML_local_improve(resultTree,model);
        logl = resultTree.compute_loglikelihood();
      }
      cout <<"Loglikelihood = " << logl << endl;
    }
    else {// BIG AML
      cout << "using heuristic for big aml"<< endl;
      seqs.clear();
      Sequence::readSequences(seqs,inSeqFile);
      //remove seqgen internal nodes 
      if ( seqgen ){
        vector<Sequence>::iterator iter = seqs.begin();
        while ( iter != seqs.end() ){
          Sequence &s = *iter;
          if ( isdigit(s.name[0]) ){
            iter = seqs.erase(iter);
          }
          else
            ++iter;
        }
      }
      Big_AML(seqs,resultTree,model);
      double logl = resultTree.compute_loglikelihood();
      cout <<"Loglikelihood = " << logl << endl;
    }
    //--------
    resultTree.recalcNodeStructure();
    resultTree.setNodeNames();    
    outTreeFile << resultTree << endl;
    outSeqFile << "   " << resultTree.getNumNodes() << "  " <<  resultTree.getRoot()->data.s.seq.length() << endl;  
    resultTree.printSequences(outSeqFile);    
    if ( verbose )
      resultTree.verbosePrint( cout );
  }
  cout << "**************"<< endl;
  
  //close out files
  outTreeFile.close();
  outSeqFile.close();

  //close in files
  inTreeFile.close();
  inSeqFile.close();
  
  CATCH_EXCEPTION();

  return 1;
}
