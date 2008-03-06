
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include "file_utils.hpp"
#include "log_utils.hpp"
#include "stl_utils.hpp"
#include <string>
#include <fstream>
#include <time.h>
#include "SequenceTree.hpp"
#include "Sequence.hpp"
#include "SequenceParsimony.hpp"
#include "LeastSquaresFit.hpp"
#include "NeighborJoining.hpp"
#include "Simulator.hpp"
#include "Sequences2DistanceMatrix.hpp"
#include "arg_utils_ext.hpp"

#include <vector>
#include "arg_utils.h"

using namespace std;

//For each model tree. we have the following trees with scores
//1. unique bootstrap trees trees, parsimony score, rfdist, l2score.
//2. Contracted trees based one parsimony trees, parsimony score, rfdist, l2score
//

typedef struct {
  int numLeafs;
  int seqLen;
  float diameterFactor;

  vector<Sequence> seqs;
  StrDblMatrix dm;
  SequenceTree speciesTree;
  SequenceTree contractedSpeciesTree;
  double speciesParsimonyScore;
  double speciesL2Score;

  //the first tree is the nj tree.
  vector<SequenceTree> bootTrees;
  vector<double> parsimonyScores;
  vector<double> rfDistances;
  vector<double> l2Scores;
  vector<double> sumEdgeLengths;
  
//   vector<SequenceTree> contractedTrees;
//   vector<double> contractedParsimonyScores;
//   vector<double> contractedRFDistances;
//   vector<double> contractedL2Scores;
  //Since edges can have length 0 the l2 score of a contracted tree
  //can not be better then the score for the uncontracted tree.
  // Therefore the l2Score is that of the best tree that was contracted into
  // contracted tree.

} TestCase;

double
computeAverage(std::vector<double> &v){
  double sum=0;
  for(size_t i=0;i<v.size();i++)
    sum+=v[i];
  return sum/v.size();
}

int
findMax(std::vector<double> &v){
  double max=-1;
  size_t maxi=99999999;
  for(size_t i=0;i<v.size();i++)
    if(v[i]>max){
      maxi=i;
      max=v[i];
    } 
  return maxi;
}
int
findMin(std::vector<double> &v){
  double min=9999999;
 int mini=-1;
  for(size_t i=0;i<v.size();i++)
    if(v[i]<min){
      mini=i;
      min=v[i];
    } 
  return mini;
}

//The minimum element in v. If two elements same size use v2 to decide min
int
findMin(std::vector<double> &v, std::vector<double> &v2){
  double min=9999999;
  int mini=-1;
  for(size_t i=0;i<v.size();i++)
    if(v[i]<=min){
      if(v[i]==min){
	if(v2[i]<v2[mini]){
	  mini=i;
	  min=v[i];
	}
      }
      else{
	mini=i;
	min=v[i];
      }
    } 
  return mini;
}
//The stats
// we compute 
//

void
generateTestCase(int numBootTrees,
		 int numLeafs,
		 int seqLen,
		 int ultrametricDeviation,
		 float diameterFactor,
		 TestCase &test){//the unique topologies){
	
  //====
  test.numLeafs = numLeafs;
  test.seqLen = seqLen;
  test.diameterFactor = diameterFactor;
  test.seqs.clear();
  test.bootTrees.clear();
  test.parsimonyScores.clear();
  test.rfDistances.clear();
  test.l2Scores.clear();
  test.sumEdgeLengths.clear();
 //  test.contractedTrees.clear();
//   test.contractedParsimonyScores.clear();
//   test.contractedRFDistances.clear();
//   test.contractedL2Scores.clear();
  //---------------------------------------
  //1.      
  //create species tree and sequences
  createRandomTreeUsingBeep(test.speciesTree,numLeafs, ultrametricDeviation);
  LINE(); PRINT(test.speciesTree); LINE();
  evolveSequencesUsingSeqGen(test.speciesTree, seqLen, diameterFactor,true);
  shortcutDegree2Nodes(test.speciesTree);

  
  //get the leaf sequences
  SequenceTree::NodeVector leafs;
  test.speciesTree.addLeafs(leafs);
  for(size_t i=0;i<leafs.size();i++){
    test.seqs.push_back(SEQUENCE(leafs[i]));
  }

  //scores of spieces tree
  test.speciesParsimonyScore = computeMostParsimoniousSequences(test.speciesTree);
  computeEdgeLengths(test.speciesTree);
  test.contractedSpeciesTree = test.speciesTree;
  int numContracted = contractEdgesShorterThan(test.contractedSpeciesTree,0);
  PRINT(numContracted);
  if(numContracted>0){
    test.speciesTree.drawTree(cout);
    test.contractedSpeciesTree.drawTree(cout);
  }


  //-------------------------
  //2.
  //skapa distance matri
  vector<DNA_b128_String> b128s;
  Sequences2DNA_b128(test.seqs,b128s);
  sequence_translation_model model;
  model.tstvratio = 2;
  model.no_tstvratio = false;
  fillMatrix_K2P(test.dm, b128s,model);
  PRINT(test.dm);
  PRINT(applyFixFactor(test.dm,2.0));

  test.speciesL2Score = computeLeastSquaresEdgeLengths(test.dm,test.speciesTree);
  
  //NJ
  SequenceTree njTree;
  computeBioNJTree(test.dm,njTree);
  PRINT(njTree);
  test.bootTrees.push_back(njTree);
  mapSequencesOntoTree(njTree,test.seqs);
  double rfDist =computeRobinsonFoulds(njTree,test.contractedSpeciesTree);  
  test.rfDistances.push_back(rfDist);
  //L 2 fitt (the contracted tree can not have a better l2 fit)
  double l2fit = computeLeastSquaresEdgeLengths(test.dm,njTree);
  test.l2Scores.push_back(l2fit);
  PRINT(l2fit);
  double sumEdges = sumOfEdgeLengths(njTree);
  PRINT(sumEdges);
  test.sumEdgeLengths.push_back(sumEdges);
  //compute parsi and collapse
  size_t parsifit = computeMostParsimoniousSequences(njTree);
  test.parsimonyScores.push_back(parsifit);
  PRINT(parsifit);


  //---------------------------
  //3, BOOT
  StrDblMatrix dm;
  for(int booti=0;booti<numBootTrees;booti++){
    LINE();PRINT(booti);
    vector<Sequence> boot;
    Sequence::bootstrapSequences(test.seqs,boot);
    b128s.clear();
    Sequences2DNA_b128(boot,b128s);
    fillMatrix_K2P(dm, b128s,model);
    PRINT(applyFixFactor(dm,2.0));
    //PRINT(dm);

    SequenceTree bootTree;
    //    computeNeighborJoiningTree(dm,bootTree);
    computeBioNJTree(dm,bootTree);
    PRINT(bootTree);

    
    double rfDist =computeRobinsonFoulds(bootTree,test.contractedSpeciesTree);  
    PRINT(rfDist);
    double l2fit = computeLeastSquaresEdgeLengths(dm,bootTree);
    PRINT(l2fit);
    computeLeastSquaresEdgeLengths(test.dm,bootTree);
    double sumEdges = sumOfEdgeLengths(bootTree);
    PRINT(sumEdges);
    
    //compute parsi and collapse
    mapSequencesOntoTree(bootTree,test.seqs);
    size_t parsifit = computeMostParsimoniousSequences(bootTree);
    PRINT(parsifit);


    bool treeExists = false;
    for(size_t i=0;i<test.bootTrees.size();i++){
      if(treesEqual(test.bootTrees[i],bootTree)){//PENDING slow fix a hashfunction for trees
	treeExists = true;
	if(test.l2Scores[i] > l2fit)
	  test.l2Scores[i] = l2fit;
	break;
      }
    }
    PRINT(treeExists);
    if( !treeExists ){//only add tree if it hasn't allready been found
      test.bootTrees.push_back(bootTree);
      test.rfDistances.push_back(rfDist);
      test.l2Scores.push_back(l2fit);
      test.parsimonyScores.push_back(parsifit);
      test.sumEdgeLengths.push_back(sumEdges);
    }
  }
  
}




typedef struct{
  int numLeafs;
  int seqLen;
  float diameterFactor;//PRINT(test.dm);  


  //the scores for the species tress
  double l2ScoreSpecies;
  double parsimonyScoreSpecies;

  double bestRF;
  double l2ScoreOfBestRF;
  double parsimonyScoreOfBestRF;

  double rfOfNJ;
  double parsimonyScoreNJ;
  double l2ScoreNJ;

  //if two solutions have the same l2 then we choose the one with best
  //parsimony.
  double bestL2;
  double rfOfBestL2;
  double parsimonyScoreBestL2;

  //if two solutions have same parsimony then we chose the one with
  //best l2
  double bestParsimony;
  double rfOfBestParsimony;
  double l2ScoreBestParsimony;

  double bestSumEdgeLengths;
  double rfOfBestSumEdgeLengths;

} Statistic;



void
computeStatsFor(size_t numDatasets,
		int numBootTrees,
		int numLeafs,
		int seqLen,
		int ultrametricDeviation,
		float diameterFactor, 
		Statistic &st){
  st.numLeafs = numLeafs;
  st.seqLen = seqLen;
  st.diameterFactor = diameterFactor;


  double sumBestRF=0;
  double sumL2ScoreBestRF=0;
  double sumParsimonyScoreBestRF=0;
  

  double sumNJRF=0;
  double sumParsimonyScoreNJ=0;
  double sumL2ScoreNJ=0;

  //species score
  double sumSpeciesL2=0;
  double sumSpeciesParsimony=0;

  double sumBestL2=0;
  double sumBestL2RF=0;
  double sumParsimonyBestL2=0;

  double sumBestParsimony=0;
  double sumBestParsimonyRF=0;
  double sumL2BestParsimony=0;

  double sumBestSumEdges=0;
  double sumBestSumEdgesRF=0;

  string baseName;
  createFileName(baseName, numLeafs, ultrametricDeviation, seqLen, diameterFactor);
  string treesFileName = "trees_"+ baseName;
  string sequencesFileName =  "sequences_"+baseName;
  string dmFileName = "dm_" +baseName;
  string bootFileName = "boottrees_" +baseName;
  ofstream treeFile;
  open_write_stream(treesFileName,treeFile);
  ofstream sequencesFile;
  open_write_stream(sequencesFileName,sequencesFile);
  ofstream dmFile;
  open_write_stream(dmFileName,dmFile);
  ofstream bootFile;
  open_write_stream(bootFileName,bootFile);

  for(size_t i=0; i<numDatasets;i++){
    TestCase test;
    generateTestCase(numBootTrees,numLeafs,seqLen,ultrametricDeviation,diameterFactor,test); 
    
    SequenceTree tmpTree(test.speciesTree);
    SequenceTree::NodeVector nodes;
    tmpTree.addInternalNodes(nodes);
    for(size_t nodei=0;nodei<nodes.size();nodei++)
      NAME(nodes[nodei])="";
    EDGE(tmpTree.getRoot())=-1;
    treeFile << tmpTree << endl;
    Sequence::printSequences(test.seqs, sequencesFile);
    dmFile << test.dm << endl;
    bootFile << "**********************" <<endl;
    bootFile << "**********************" <<endl;
    bootFile << "**********************" <<endl;
    for(size_t booti=0;booti<test.bootTrees.size();booti++)
      bootFile << test.bootTrees[booti] << endl;
    

    int tmpi = findMin(test.rfDistances);
    sumBestRF+=test.rfDistances[tmpi];
    sumL2ScoreBestRF += test.l2Scores[tmpi];
    sumParsimonyScoreBestRF += test.parsimonyScores[tmpi];


    //the first tree is the NJ tree.
    sumNJRF += test.rfDistances[0];
    sumParsimonyScoreNJ += test.parsimonyScores[0];
    sumL2ScoreNJ += test.l2Scores[0];

    //species
    sumSpeciesL2 += test.speciesL2Score;
    sumSpeciesParsimony += test.speciesParsimonyScore;

    //RF of BEST L2
    tmpi = findMin(test.l2Scores,test.parsimonyScores);
    sumBestL2 += test.l2Scores[tmpi];
    sumBestL2RF += test.rfDistances[tmpi];
    sumParsimonyBestL2 += test.parsimonyScores[tmpi];
    
    //RF of BEST PARSIMONY
    tmpi = findMin(test.parsimonyScores,test.l2Scores);
    sumBestParsimony += test.parsimonyScores[tmpi];
    sumBestParsimonyRF += test.rfDistances[tmpi];
    sumL2BestParsimony += test.l2Scores[tmpi];
   
    //RF of BEST EDGE LENGTHS
    tmpi = findMin(test.sumEdgeLengths,test.l2Scores);
    sumBestSumEdges += test.sumEdgeLengths[tmpi];
    sumBestSumEdgesRF += test.rfDistances[tmpi];
  }
  
  treeFile.close();
  sequencesFile.close();
  dmFile.close();
  bootFile.close();

  //set the average value over all the runs
  st.bestRF = sumBestRF/numDatasets;
  st.l2ScoreOfBestRF = sumL2ScoreBestRF/numDatasets;
  st.parsimonyScoreOfBestRF = sumParsimonyScoreBestRF/numDatasets;

  st.rfOfNJ = sumNJRF/numDatasets;
  st.parsimonyScoreNJ = sumParsimonyScoreNJ/numDatasets;
  st.l2ScoreNJ = sumL2ScoreNJ/numDatasets;

  st.l2ScoreSpecies = sumSpeciesL2/numDatasets;
  st.parsimonyScoreSpecies = sumSpeciesParsimony/numDatasets;

  st.bestL2 = sumBestL2/numDatasets;
  st.rfOfBestL2 = sumBestL2RF/numDatasets;
  st.parsimonyScoreBestL2 = sumParsimonyBestL2/numDatasets;

  st.bestParsimony = sumBestParsimony/numDatasets;
  st.rfOfBestParsimony = sumBestParsimonyRF/numDatasets;
  st.l2ScoreBestParsimony = sumL2BestParsimony/numDatasets;

  st.bestSumEdgeLengths = sumBestSumEdges/numDatasets;
  st.rfOfBestSumEdgeLengths = sumBestSumEdgesRF/numDatasets;
}




//WE CREATE FOUR DIFFERENT PLOTS FOR EACH FIXED TREE SIZE

//1. the rf distance for the different sequence lengths for each
//optimal tree of each criteria.

//2. the parsimony score of the best from each criteria.

//3. the l2 score of the best of each criteria.

//4. The different scores of the best tree.
void
createPlotsForTreeSize(int numLeafs, double diameterFactor, 
		       std::vector<Statistic> &stats){
  
  //open the files
  ofstream rfFile;
  char tmpStr[1000];
  sprintf(tmpStr,"RFDistance_leafs_%d_diamF_%f.dat", numLeafs,diameterFactor);
  string rfFileName(tmpStr);
  open_write_stream(tmpStr,rfFile);
  rfFile << "#RF distance to contracted species tree\n"
	 << "#TREES WITH LEAFS = " << numLeafs<<endl
	 << "#diameterFactor = " << diameterFactor << endl;
  rfFile << "SeqLen \tBest  \tNJ  \tL2  \tParsi  \tMinEvol" << endl;

  ofstream parsimonyFile;
  sprintf(tmpStr,"ParsimonyScores_leafs_%d_diamF_%f.dat", numLeafs,diameterFactor);
  string parsimonyFileName(tmpStr);
  open_write_stream(tmpStr,parsimonyFile);
  parsimonyFile << "#Parsimony Scores for different trees\n"
		<< "#TREES WITH LEAFS = " << numLeafs<<endl
		  << "#diameterFactor = " << diameterFactor << endl;
  parsimonyFile << "SeqLen \tSpecies \tBioNJ \tL2 \tParsi" << endl;
  
  ofstream l2File;
  sprintf(tmpStr,"L2Scores_leafs_%d_diamF_%f.dat", numLeafs,diameterFactor);
  string l2FileName(tmpStr);
  open_write_stream(tmpStr,l2File);
  l2File << "#L2 Scores for different trees\n"
		  << "#TREES WITH LEAFS = " << numLeafs<<endl
	 << "#diameterFactor = " << diameterFactor << endl;
  l2File << "SeqLen \tSpecies \tBioNJ \tL2 \tParsi" << endl;
  
  ofstream bestFile;
  sprintf(tmpStr,"BestTree_leafs_%d_diamF_%f.dat", numLeafs,diameterFactor);
  string bestFileName(tmpStr);
  open_write_stream(tmpStr,bestFile);
  bestFile << "#The different scores for the tree with best RF distance\n"
		  << "#TREES WITH LEAFS = " << numLeafs<<endl
	   << "#diameterFactor = " << diameterFactor << endl;
  bestFile << "SeqLen \tRF \tL2Score \tParsi" << endl;
  
  
  //WRITE THE FILES
  for(size_t i=0;i<stats.size();i++){
    Statistic &st = stats[i];
    if(st.numLeafs!= numLeafs || st.diameterFactor!= diameterFactor)
      continue;
    
    rfFile << st.seqLen << " \t"<< st.bestRF <<  " \t" <<st.rfOfNJ << 
      " \t"<<st.rfOfBestL2 << " \t" << st.rfOfBestParsimony << " \t"<< st.rfOfBestSumEdgeLengths << endl;

    parsimonyFile << st.seqLen << " \t"<< st.parsimonyScoreSpecies << " \t"<< st.parsimonyScoreNJ << " \t"
		  << st.parsimonyScoreBestL2 << " \t" << st.bestParsimony << endl;

    l2File << st.seqLen << " \t"<< st.l2ScoreSpecies << " \t"<< st.l2ScoreNJ << " \t"<< st.bestL2<< " \t" 
		  << st.l2ScoreBestParsimony << endl;

    bestFile<< st.seqLen << " \t"<<st.bestRF << " \t" << st.l2ScoreOfBestRF << " \t" << st.parsimonyScoreOfBestRF << endl;   
  }

  //close the files
  rfFile.close();
  parsimonyFile.close();
  l2File.close();
  bestFile.close();

  //CREATE THE GNUPLOTS for the three dat files.

  // rfFile << "SeqLen \tNJ \tNJCont. \tL2 \tL2Cont. \tParsi \tParsiCont." << endl;
  vector<int> columns;
  vector<string> names;
  columns.push_back(2);names.push_back("Best");
  columns.push_back(3);names.push_back("BioNJ");
  columns.push_back(4);names.push_back("L2");
  columns.push_back(5);names.push_back("Parsimony");
  columns.push_back(6);names.push_back("MinEvol");
  char title[1000];
  sprintf(title,"RF for leafs=%d and diamF=%f",numLeafs,diameterFactor);
  createGnuplotFromDatFile(rfFileName.c_str(), "SeqLen", "RF", columns,names,
			   title);

  //parsimonyFile << "SeqLen \tSpecies \tNJ \tL2 \tParsi" << endl;
  columns.clear();names.clear();
  columns.push_back(2);names.push_back("Species");
  columns.push_back(3);names.push_back("BioNJ");
  columns.push_back(4);names.push_back("L2");
  columns.push_back(5);names.push_back("Parsimony");

  sprintf(title,"Parsimony for leafs=%d and diamF=%f",numLeafs,diameterFactor);
  createGnuplotFromDatFile(parsimonyFileName.c_str(), "SeqLen", "Parsimony", columns,names,
			   title);
  
  //l2File << "SeqLen \tSpecies \tNJ \tL2 \tParsi" << endl;
  columns.clear();names.clear();
  columns.push_back(2);names.push_back("Species");
  columns.push_back(3);names.push_back("BioNJ");
  columns.push_back(4);names.push_back("L2");
  columns.push_back(5);names.push_back("Parsimony");
  sprintf(title,"L2 for leafs=%d and diamF=%f",numLeafs,diameterFactor);
  createGnuplotFromDatFile(l2FileName.c_str(), "SeqLen", "L2", columns,names,
			   title);
  
}



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
    " --numboots int           100    The number of bootstraps for each dataset\n"
    
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
  

  int numBoots = 100;
  tmp = GET_OPTION_VAL("--numboots");
  if( tmp!=NULL )
    numBoots = atoi(tmp);
  

  //CREATE THE STATISTICS
  int totalNumberOfCases = leafSizes.size()*seqLens.size()*diameterFactors.size();
  int currentCase =0;
  vector<Statistic> stats;

  for(size_t leafi=0 ; leafi<leafSizes.size(); leafi++){
    for( size_t seqi=0 ; seqi<seqLens.size() ; seqi++){
      for( size_t diami=0 ; diami<diameterFactors.size() ; diami++){
	currentCase++;
	SEPARATOR();SEPARATOR();PRINT(currentCase); PRINT(totalNumberOfCases);
	Statistic stat;
	stats.push_back(stat);
	
	computeStatsFor(numTrees,numBoots,leafSizes[leafi], seqLens[seqi],4,diameterFactors[diami],
			stats[stats.size()-1]);
      }
    } 
  }

  //create plots
  for(size_t leafi=0 ; leafi<leafSizes.size(); leafi++){
      for( size_t diami=0 ; diami<diameterFactors.size() ; diami++){
	createPlotsForTreeSize(leafSizes[leafi], diameterFactors[diami], stats);
      }
  } 

  
  
  return 1;
}












