#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>

using namespace std;

#define STRLEN 100000
#define R8S "~/10giga_volume/r8s/r8s1.71/src/r8s"
#define R8S_IN "r8s_infile.txt"
#define R8S_OUT "r8s_outfile.txt"
#define SEQGEN "~/10giga_volume/Seq-Gen.v1.3.2/source/seq-gen"
#define SEQGEN_IN "seqgen_infile.txt"
#define SEQGEN_OUT "seqgen_outfile.txt"
#define TESTDIR "testfiles/"

float
rand_float( float a, float b) {
  return a + (((float)rand())/RAND_MAX) * ( b - a );
}

float mutationProb;

void
createTest(int numLeafs,
           int seqLen,
           char *seqfile_results,
           char *treefile_results){
  
  //*********************
  //CREATE r8s in file
  ofstream r8s_infile;
  remove(R8S_IN);
  remove(R8S_OUT);
  
  r8s_infile.open(R8S_IN);
  if (!r8s_infile) {
    cout << "Could not open.\n";
    exit(1);
  }

  r8s_infile << "#nexus" << endl;
  r8s_infile << "begin trees;" << endl;
  r8s_infile << "end;" << endl;
  r8s_infile << "begin r8s;" << endl;
  r8s_infile << "simulate diversemodel=bdback seed="<< rand() << " ntaxa=" << numLeafs << 
    "maxrate=0.20" //isaac added to test for short edges
             <<" T=0 ;" << endl;
  r8s_infile << "describe tree=0 plot=tree_description;" << endl;
  r8s_infile << "end;" << endl;
  r8s_infile.close();

  //execute r8s
  cout<< " about to execute: " << R8S " -b -f " R8S_IN " > " R8S_OUT << endl;
  system(R8S " -b -f " R8S_IN " > " R8S_OUT );  
  //*********************
  //Read the tree string
  ifstream r8s_outfile;
  r8s_outfile.open(R8S_OUT);
  if (!r8s_outfile) {
    cout << "Input file cannot be opened.\n";
    exit(1);
  }
  char str[STRLEN];
  
  bool SIMTREEfound = false;
  while ( true) {
    r8s_outfile >> str;
    if ( strcmp(str,"SIMTREE") == 0 )
      SIMTREEfound = true;
    
    if ( SIMTREEfound ){
      if ( strcmp(str,"=") == 0 ){
        r8s_outfile >> str;
        break; // the TREE IS NOW IN STR
      }
    }
  }
  r8s_outfile.close();
  //cout << "the ultrametric tree is: " << str << endl;
  //****************
  //CHANGE THE EDGE LENGTHS with ultrametricDeviation
  char tree[STRLEN];
  char *t = tree;
  char *s = str;
  while ( *s != '\0' ){
    //cout << "cpying: " << *s << endl;
    if ( *s == ':' ){
      *t = *s;
      t++; s++;      
      float elen = 0;
      sscanf(s,"%f",&elen);
      //cout << "float is: " << elen << endl;
      while ( isdigit(*s) || *s == '.' )
        s++;

      //old code used for FNJ
      //float new_elen = elen*rand_float(((float)1)/ultrametricDeviation, 
      //       ultrametricDeviation);

      //the length is a random number in the interval 0-0.25
      float new_elen = rand_float(0,mutationProb);
      
      int numWritten = sprintf(t,"%f",new_elen);
      t = t+numWritten;
    }
    else {
      *t = *s;
      t++; s++;
    }
  }
  *t = '\0';

  //cout << "THE TREE IS: " << tree << endl;
  //*******************************
  //Write to seq-gen file exucute SEQ-GEN
  ofstream seqgen_infile;
  remove(SEQGEN_IN);
  remove(SEQGEN_OUT);
  
  seqgen_infile.open(SEQGEN_IN);
  if (!seqgen_infile) {
    cout << "Could not open.\n";
    exit(1);
  }

  seqgen_infile << tree << endl;
  seqgen_infile.close();

  //execute seq-gen
  sprintf(str,  SEQGEN " -l%d "
                    " -wa " //write ancestor sequences  
          //JUKES CANTOR
          //	  " -s %f -n 1 -m HKY -t 0.5 -f  0.25 0.25 0.25 0.25 < "
          " -n1 -mHKY -t0.5 -f0.25,0.25,0.25,0.25 < " //without diamter factor
          //K2P with fix ratio 2
          //" -s %f -n 1 -m HKY -t 2 -f  0.25 0.25 0.25 0.25 < "
          SEQGEN_IN " > " SEQGEN_OUT,
          seqLen);//, diameterFactor);
  cout << "The command line is: " << str << endl;
  system(str);

  //********************
  //Create test file
  //write the model tree
  sprintf(str, "cat " SEQGEN_IN " >> %s", treefile_results);
  system(str);
  //write the sequences
  sprintf(str, "cat " SEQGEN_OUT " >> %s", seqfile_results);
  system(str);
}


int
main(int argc, char **argv){
  cout << "Run as ./generate_phylogenies 0.25" << endl;
  mutationProb = atof(argv[1]);
  int seqLengths[] = {10,50,100,500,1000,6000};
  int numSL = 6;
  int leafSizes[] = {5,10,15,20,30,40,50,60,70,80,90,100};
  int numLS = 12;
  int numTreesPerSize = 20;

  char treefilename[100];
  char seqfilename[100];
  //remove all test files
  system("rm -f " TESTDIR "test_leafs*" );

  int test = 1;
  //leaf sizes
  for ( int ls = 0 ; ls < numLS ; ls++ ){
    for ( int sl = 0 ; sl < numSL ; sl++ ){
      //num trees for each combination
      sprintf(treefilename,TESTDIR "trees_leafs_%d_sqlen_%d_mprob_%f.txt",
              leafSizes[ls],seqLengths[sl],mutationProb);
      sprintf(seqfilename,TESTDIR "sequences_leafs_%d_sqlen_%d_mprob_%f.txt",
              leafSizes[ls],seqLengths[sl],mutationProb);
      //char numtests[100];
      //sprintf(numtests,"echo %d    num trees in file > %s", numTreesPerSize,filename);
      //system(numtests);
      for ( int tps = 0 ; tps < numTreesPerSize ; tps++ ){
        cout <<"*************************** Doing test " << test << endl;
        test++;
	
        createTest(leafSizes[ls], seqLengths[sl], seqfilename,treefilename);
      }
    }
  }
}






































