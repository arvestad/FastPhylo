
#include <math.h>
#include "arg_utils_ext.hpp"
#include "stl_utils.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include "SequenceTree.hpp"
#include <algorithm>  // Include algorithms
#include <functional>
#include "stl_utils.hpp"
#include <list>

using namespace std;
typedef __gnu_cxx::hash_set<SequenceTree::Node*, objhash_ptr, objeq_ptr> node_set;
typedef __gnu_cxx::hash_map<Sequence, std::vector<Sequence>,objhash, objeq> SequenceVector_map;

#define STRLEN 100000
#define BEEP "~/10giga_volume/beep/beep_generateTree"
#define BEEP_OUT "beep_outfile.txt"

#define ROSE "~/10giga_volume/rose-1.3/src/rose"
#define ROSE_OUT "rose_outfile.txt"
#define ROSE_IN "rose_infile.txt"

#define TESTDIR "testfiles/"

string
gaplessString(string &seq);
int
chooseRandomPosition(list<int> &positions);
void
removeAllIntsInInterval(list<int> &positions, int start, int end);

float
rand_float( float a, float b) {
  return a + (((float)rand())/RAND_MAX) * ( b - a );
}
int
rand_int( int a, int b){
  return a+(int) (b*1.0*rand()/(RAND_MAX+1.0));
}

list<char>::iterator
getIteratorAtPos(list<char> &l, 
		 size_t pos){
  list<char>::iterator iter = l.begin();
  for( size_t i=0 ; i<pos ; i++ )
    ++iter; 
  return iter;
}

struct less_root_dist : public binary_function<SequenceTree::Node *, SequenceTree::Node *, bool> {
  bool operator()(SequenceTree::Node *x, SequenceTree::Node *y) { return x->data.dbl < y->data.dbl; }
};



void
generateSequencesUsingRose(SequenceTree &tree,int seqlen, double exp_sub, double exp_indel){
  //ROSE
  //generate sequences according to base tree
  ofstream rose_parameter_file;
  open_write_stream(ROSE_IN,rose_parameter_file);

  //the root to leaf distance is 100
  double subrate = (-log(((1.0-exp_sub)-0.25)/0.75))/100.0;
  float threshold = exp_indel/100.0;
  
  rose_parameter_file << 
    "# rose parameters\n"
    "StdOut = False\n"
    "OutputFilebase = \"" ROSE_OUT "\"\n"
    "InputType = 4\n" // DNA
    "TheAlphabet = \"ACGT\"\n"
    "TheFreq = [.25,.25,.25,.25]\n"
    
    "TheInsertThreshold = " << fixed << threshold << "\n"
    "TheDeleteThreshold = " << fixed << threshold << "\n"
    "TheInsFunc = [.2,.2,.2,.1,.1,.1,.1]\n"
    "TheDelFunc = [.2,.2,.2,.1,.1,.1,.1]\n"
    
    "TheDNAmodel = \"JC\"\n"
    //"MeanSubstitution = 0.01342302\n" // 1 percent mutations for t=1
    //the distance in the tree from root to leaf is 100.
    "MeanSubstitution = " << subrate << "\n" 
    "TransitionBias = 1.0\n" // reduces K2P to JC
    "TTratio = 0.0\n" // reduces F84 to JC
    "SequenceNum = " << tree.getNumLeafs() << endl <<//all leafs i.e. numLeafs of tree
    "ChooseFromLeaves = True\n" //output leaf sequences
    //root sequence TheSequence = "AGTCTGTACTATAATGGGAGGAAAGCC"
    "SequenceLen = " << seqlen << endl<<
    "TheTree = " << tree << endl
    ;

  rose_parameter_file.close();

  string rose_str = ROSE;  
  rose_str += " " ROSE_IN;

  cout<< " about to execute: " << rose_str<< endl;
  system(rose_str.c_str() );  

  //---
  //map sequences onto tree
  ifstream seqin;
  open_read_stream(ROSE_OUT ".phy", seqin);
  mapSequencesOntoTree(tree,seqin);
  seqin.close();

}


void
creatRandomTreeUsingBeep(SequenceTree &outTree, int numLeafs){
  //*********************
  //RUN beep in file
  remove(BEEP_OUT);
  
  //execute beep
  //creates an ultrametric birth death tree with hight 1.0.
  string beep_str = BEEP;
  beep_str += string(" -o " BEEP_OUT " -To params.true -Gn -Gt -Bp 0.4 0.4 ");
  beep_str = beep_str+numLeafs;
  
  cout<< " about to execute: " << beep_str<< endl;
  system(beep_str.c_str() );  
  //*********************
  //Read the tree string
  ifstream beep_outfile;
  beep_outfile.open(BEEP_OUT);
  if (!beep_outfile) {
    cout << "Input file cannot be opened.\n";
    exit(1);
  }
  char str[STRLEN];
  beep_outfile >> str;//the tree is now in str    
  beep_outfile.close();


  //base tree from beep the root to leaf path has length 1.

  SequenceTree tmpTree(str);
  tmpTree.getRoot()->data.dbl = -1;
  
  SequenceTree::NodeVector nodes;
  tmpTree.addNodesInPrefixOrder(nodes);
  //scale edges
  for(size_t i=1 ; i<nodes.size() ; i++){//skip root
    nodes[i]->data.dbl = nodes[i]->data.dbl*100; //the distance from root2leaf is 100
  }
  
  outTree = tmpTree;
}

void
performOneTransferOnTree(SequenceTree &tree){

  SequenceTree::NodeVector nodes;
  
  //change so that each node has its distance to the root
  nodes.clear();
  tree.addNodesInPrefixOrder(nodes);
  nodes[0]->data.dbl = 0;//root has distance 0 to itself
  for(size_t i=1 ; i<nodes.size() ; i++ ){
    nodes[i]->data.dbl += nodes[i]->getParent()->data.dbl;
  }

  //sort nodes according to their distance to the root
  sort(nodes.begin(), nodes.end(),less_root_dist());

  //find all concurrent edges (transfers can only be happen between
  //edges that are alive at the same time)
  vector<SequenceTree::NodeVector> concurrentEdges; 
  vector<double> evolutionaryTimeIncommon;  
  ASSERT_EQ(nodes[0], tree.getRoot());
  
  SequenceTree::Node *lastNode = tree.getRoot();
  node_set currentEdges;
  currentEdges.insert(tree.getRoot()->getRightMostChild());
  currentEdges.insert(tree.getRoot()->getLeftMostChild());
  double totalEvolutionaryTime=0;
  for(size_t i=1 ; i<nodes.size() ; i++){
    double timeInCommon = currentEdges.size()*(nodes[i]->data.dbl-lastNode->data.dbl);
    totalEvolutionaryTime += timeInCommon;
    //PRINT(timeInCommon);
    evolutionaryTimeIncommon.push_back(timeInCommon);
    SequenceTree::NodeVector edges(currentEdges.begin(),currentEdges.end());
    concurrentEdges.push_back(edges);
    
    if(nodes[i]->isLeaf() )
      break;

    lastNode = nodes[i];
    currentEdges.erase(nodes[i]);
    currentEdges.insert(nodes[i]->getRightMostChild());    
    currentEdges.insert(nodes[i]->getLeftMostChild());
  }

  //FIND the transfer:
  //the transfers go from the edge above startEdge to the edge above
  //endEdges and is on distance timeOfTransfer from the root.

 new_rand_time:
  float randTotalEvolTimeBeforeTransfer = rand_float(0,totalEvolutionaryTime);
  if( randTotalEvolTimeBeforeTransfer<=evolutionaryTimeIncommon[0] ){
    goto new_rand_time;//cannot have transfer between the two children of the root
  }
  double timeSum=0;
  for(size_t j=0 ; j<evolutionaryTimeIncommon.size() ;j++ ){
    if( randTotalEvolTimeBeforeTransfer>timeSum+evolutionaryTimeIncommon[j] ){
      timeSum+=evolutionaryTimeIncommon[j];
      continue;
    }
      
    //randTotalEvolTimeBeforeTransfer<=tmpTime+evolutionaryTimeIncommon[j]
    //tansfer position found
    
    //choose the two edges uniformly from the concurrentEdges
    int numConcurrent = concurrentEdges[j].size();
    int startEdge = rand_int(0,numConcurrent-1); 
    if(numConcurrent==3){//special fix for the case when there are only three concurrent
      for(int i=0 ; i<numConcurrent ; i++)
	if(concurrentEdges[j][i]->getParent()->isRoot() )
	  startEdge = i;
    }
    int endEdge; //edge with different parent and not the root as parent
    do{
      endEdge = rand_int(0,numConcurrent-1);
    }while(concurrentEdges[j][startEdge]->getParent()==concurrentEdges[j][endEdge]->getParent() ||
	   concurrentEdges[j][endEdge]->getParent()->isRoot());
      
    //choose the time of the event uniformly from the common interval
    double commonIntervalLength = evolutionaryTimeIncommon[j]/concurrentEdges[j].size();
    double timeOfEvent = rand_float(0,commonIntervalLength);
    timeOfEvent += nodes[j]->data.dbl;
    
    SequenceTree::Node *start = concurrentEdges[j][startEdge];
    SequenceTree::Node *end = concurrentEdges[j][endEdge];
    SEPARATOR();
    PRINT(timeOfEvent);
    PRINT_EXP(tree.drawSubtree(cout,end));
    PRINT("to"); 
    PRINT_EXP(tree.drawSubtree(cout,start));
    SEPARATOR();

    SequenceTree::Data_type data;
    data.dbl = timeOfEvent;
    tree.insertNodeOnPathToParent(start,data);
    SequenceTree::Node *oldParent = end->getParent();
    tree.moveNode(end, start->getParent());
    tree.collapse(oldParent);
    break; //the transfer has been found and performed
  }//find one transfer position loop

  
  //Change the edge lengths so that they are the regular lengths
  //and not the node to root distance.
  nodes.clear();
  tree.addNodesInPrefixOrder(nodes);
  for(size_t i=nodes.size()-1 ; i!=0 ; i-- )
    nodes[i]->data.dbl -= nodes[i]->getParent()->data.dbl;
  nodes[0]->data.dbl = -1;
}

void
createTest(string testname,
	   int numLeafs, 
	   double exp_sub, 
	   double exp_indel, 
	   int numTransferTrees, 
	   int numTransfersPerTree,
	   int seqLenSpeciesTree,
	   int seqLenTransferTrees){
  
  //-------------------------------
  // Create the species tree
  SequenceTree speciesTree;
  creatRandomTreeUsingBeep(speciesTree,numLeafs);

  SEPARATOR();PRINT(speciesTree);speciesTree.drawTree(cout);

  //----------------------------------------
  //create the transfer trees
  std::vector<SequenceTree> transferTrees;
  for(int i=0 ; i<numTransferTrees ; i++ ){
    SequenceTree transferTree(speciesTree);
    for(int j=0 ; j<numTransfersPerTree ; j++){
      performOneTransferOnTree(transferTree);
      SEPARATOR();PRINT(transferTree);transferTree.drawTree(cout);
    } 
    transferTrees.push_back(transferTree);
    SEPARATOR();PRINT(transferTree);transferTree.drawTree(cout);
  }  
  
  SequenceTree::NodeVector nodes;
  speciesTree.addLeafs(nodes);

  //a map from the leaf name to all sequences in that leaf.
  SequenceVector_map sequencesInLeaf;
  for(size_t i=0 ; i<nodes.size() ; i++ ){
    vector<Sequence> seqVec;
    seqVec.push_back(nodes[i]->data.s);
    sequencesInLeaf[(nodes[i]->data.s)]=seqVec;
  }

  //----------------------------------
  //create the sequences and write to files
  
  //species tree
  generateSequencesUsingRose(speciesTree,seqLenSpeciesTree,exp_sub,exp_indel);
  ofstream speciesTreeFile;
  string speciesFileName = testname +"_species";
  //tree file
  open_write_stream((speciesFileName+".tree").c_str(), speciesTreeFile);
  speciesTreeFile << speciesTree << endl;
  speciesTreeFile.close();
  //alignment file
  ofstream speciesAlignFile;
  open_write_stream((speciesFileName+".align").c_str(), speciesAlignFile);
 
  printSequencesPhylip(nodes,speciesAlignFile);
  speciesAlignFile.close();
  
  //enter the generated sequences
  for(size_t i=0 ; i<nodes.size() ; i++ ){
    SequenceVector_map::iterator iter=sequencesInLeaf.find(nodes[i]->data.s);
    if( iter==sequencesInLeaf.end() ){
      PROG_ERROR("couldn't find the leaf...it should be there");
    }
    (*iter).second.push_back(nodes[i]->data.s);
  }


  //transfer trees
  for(int i=0 ; i<numTransferTrees ; i++ ){
    SequenceTree &transferTree = transferTrees[i];
    generateSequencesUsingRose(transferTree,seqLenTransferTrees,exp_sub,exp_indel);
    string transferFileName = testname +"_transfer";
    transferFileName = transferFileName+i;
    //tree file
    ofstream out;
    open_write_stream((transferFileName+".tree").c_str(), out);
    out << transferTree << endl;
    out.close();
    //alignment file
    open_write_stream((transferFileName+".align").c_str(), out);
    nodes.clear();
    transferTree.addLeafs(nodes);
    printSequencesPhylip(nodes,out);
    out.close();

    //enter the generated sequences
    for(size_t i=0 ; i<nodes.size() ; i++ ){
      SequenceVector_map::iterator iter=sequencesInLeaf.find(nodes[i]->data.s);
      if( iter==sequencesInLeaf.end() ){
	PROG_ERROR("couldn't find leaf...it should be there");
      }
      (*iter).second.push_back(nodes[i]->data.s);
    }
  }

  //combined alignment file  
  ofstream combinedAlignFile;
  open_write_stream((testname+"_combined.align").c_str(), combinedAlignFile);
  vector<Sequence>::iterator veciter = (*sequencesInLeaf.begin()).second.begin();
  int alignlen=0;
  for(; veciter !=(*sequencesInLeaf.begin()).second.end() ; ++veciter ){
    alignlen += (*veciter).seq.length(); 
  }
  combinedAlignFile << numLeafs << "\t " << alignlen << endl;
  SequenceVector_map::iterator iter = sequencesInLeaf.begin();
  for( ; iter != sequencesInLeaf.end() ; ++iter){
    veciter = (*iter).second.begin();
    combinedAlignFile << (*veciter);//output leafName and sequence of species tree
    ++veciter;
    for( ; veciter!=(*iter).second.end() ; ++veciter ){//output the sequences of the tranfer blocks
      combinedAlignFile << "  " << (*veciter).seq;
    }
    combinedAlignFile << endl;
  }
  combinedAlignFile.close();


  //random combination file
  ofstream combinedRandFile;
  open_write_stream((testname+"_combinedRand.seqs").c_str(), combinedRandFile);  
  ofstream combinedIndexFile;
  open_write_stream((testname+"_combinedRand.index").c_str(), combinedIndexFile);  
  
  iter = sequencesInLeaf.begin();
  for( ; iter!= sequencesInLeaf.end() ; ++iter){
    veciter = (*iter).second.begin();//the first is only the name
    combinedRandFile << ">" << (*veciter).name << endl;
    combinedIndexFile << ">" << (*veciter).name << endl;
    //a list of all the characters that makes up the sequences
    list<char> seqlist;
    //add the species sequence to the list
    list<char>::iterator seqiter = seqlist.end();
    ++veciter;//the second is the species sequence
    PRINT((*veciter).seq);PRINT((*veciter).name);
    string gapless = gaplessString((*veciter).seq);
    PRINT(gapless);
    seqlist.insert(seqiter, gapless.begin(), gapless.end());
    PRINT(gapless.size());

    list<int> availablePositions;
    for(int pos=seqLenTransferTrees+1; pos<((int)gapless.length())-seqLenTransferTrees-1 ; pos++ )
      availablePositions.push_back(pos);
    LINE();PRINT(availablePositions.size());

    //find random positions for the other sequences to be inserted in
    vector<list<char>::iterator> randomPositions;
    vector<int> randomInts;

    ++veciter;//skip sequence that evolved according to species tree
    for( ; veciter!=(*iter).second.end() ; ++veciter ){      
      int randpos = chooseRandomPosition(availablePositions);
      removeAllIntsInInterval(availablePositions, randpos-seqLenTransferTrees,randpos+seqLenTransferTrees);
      LINE();PRINT(availablePositions.size());
      randomInts.push_back(randpos);
      randomPositions.push_back(getIteratorAtPos(seqlist,randpos));
    }

    //insert the transfer blocks  in their positions
    
    vector<list<char>::iterator>::iterator positer = randomPositions.begin();
    vector<int>::iterator intiter = randomInts.begin();
    veciter = (*iter).second.begin();
    ++veciter; ++veciter;//skip the empty seq and the species seq
    for( ; veciter!=(*iter).second.end() ; ++veciter,++positer,++intiter ){
      gapless = gaplessString((*veciter).seq);
      PRINT(gapless);
      seqlist.insert(*positer,gapless.begin(), gapless.end());

      //update positions
      int position = *intiter;  
      for(size_t i=0 ; i<randomInts.size() ; i++ ){
	if(randomInts[i]>position){
	  randomInts[i] += gapless.size();
	}
      }
    }
    assert(positer==randomPositions.end());

    //output the sequence without gaps
    int outputedChars=0;
    list<char>::iterator chariter =seqlist.begin();
    for( ; chariter!=seqlist.end() ; ++chariter ){
      char c = *chariter;
   
      outputedChars++;
      combinedRandFile << c;
      if( outputedChars%50 == 0 ){
	combinedRandFile << endl;
      }
    }
    combinedRandFile << endl << endl;
    //output the indexes of the blocks
    veciter = (*iter).second.begin();
    ++veciter; ++veciter;//skip species seq
    for(size_t i=0;i<randomInts.size();i++,++veciter){
      gapless = gaplessString((*veciter).seq);
      combinedIndexFile << "block_"<<i<<"  [" <<randomInts[i] <<"," << (randomInts[i]+gapless.length()) << "]"<<endl;
    }
  }
  combinedRandFile.close();
  combinedIndexFile.close();
}

//returns the number of characters that were inserted
string
gaplessString(string &seq){

  string tmp;
  for(size_t i=0 ; i<seq.length() ; i++){
    if(seq[i]!='-'){
      tmp.append(1,seq[i]);
    }
  }

  return tmp;
}

int
chooseRandomPosition(list<int> &positions){
  int randpos = rand_int(0,positions.size());
  if(positions.size()==0){
    USER_ERROR("not enough positions");
  }
  list<int>::iterator iter = positions.begin();
  while(randpos>0){
    ++iter;
    randpos--;
  }
  return *iter;
}

void
removeAllIntsInInterval(list<int> &positions, int start, int end){
  list<int>::iterator iter = positions.begin();
  while( iter!=positions.end() ){
    if(*iter>start && *iter<end ){
      iter = positions.erase(iter);
    }
    else{
      ++iter;
    }
  }
}



void
print_options(char *note = NULL){
  if ( note != NULL ){
    cout << "ERROR: " << note <<endl<< endl;
  }
  cout <<
    "OPTIONS              MAN/DEF   DESCRIPTION\n"
    " -fileprefix str     dataset   The prefix of all output files\n"
    " -num_datasets int       1     The number of datasets to generate.\n"
    " -num_leafs int         10     The number of leafs the trees should have.\n"
    " -species_seq_len int  100     The length of the sequences evolved according to species tree.\n"
    " -transfer_trees int     3     Number of blocks that have evolved according to other tree.\n"
    " -transfers_per_tree     2     The number of transfers in the species tree for each block\n"
    " -block_len int         20     The length of the blocks.\n"
    " -exp_sub flt          0.2     The expected number substitutions in one position from root the root down to leafs\n"
    " -exp_indel flt       0.01     The expected number indels in one position from root the root down to leafs\n"
    "\n"
    " -h or --help               Print this help message\n";
  cout << endl;

  exit(1);
}




int
main(int argc, char **argv){


  if ( HAS_OPTION("-h") || HAS_OPTION("--help") )
    print_options();

  
  int num_datasets=1;
  SET_INT_OPTION_VAL("-num_datasets",&num_datasets);

  string fileprefix("dataset");
  char *prefix= GET_OPTION_VAL("-fileprefix");
  if(prefix != NULL ) 
    fileprefix = string(prefix);

  int num_leafs=10;
  SET_INT_OPTION_VAL("-num_leafs",&num_leafs);

  float exp_sub=0.2; 
  SET_FLOAT_OPTION_VAL("-exp_sub",&exp_sub);
  if(exp_sub > 0.75 )
    print_options("expected number of changes can not be more than 0.75");
  float exp_indel=0.01; 
  SET_FLOAT_OPTION_VAL("-exp_indel",&exp_indel);
  

  int numTransferTrees=3;
  SET_INT_OPTION_VAL("-transfer_trees",&numTransferTrees);

  int numTransfersPerTree=2;
  SET_INT_OPTION_VAL("-transfers_per_tree", &numTransfersPerTree);
  
  int seqLenSpeciesTree=100;
  SET_INT_OPTION_VAL("-species_seq_len", &seqLenSpeciesTree);

  int seqLenTransferTrees=20;
  SET_INT_OPTION_VAL("-block_len", &seqLenTransferTrees);


  for(int i=1 ; i<=num_datasets ; i++){
    createTest((fileprefix+"_")+i,num_leafs,exp_sub,exp_indel,numTransferTrees,numTransfersPerTree,
	       seqLenSpeciesTree,seqLenTransferTrees);
  }
}







