#include "stl_utils.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>

using namespace std;


#define STRLEN 100000
#define BEEP "~/10giga_volume/beep/beep_generateTree"
#define BEEP_OUT "beep_outfile.txt"
#define SEQGEN "~/10giga_volume/Seq-Gen.v1.3.2/source/seq-gen"
#define SEQGEN_IN "seqgen_infile.txt"
#define SEQGEN_OUT "seqgen_outfile.txt"
#define TESTDIR "testfiles/"

float
rand_float( float a, float b) {
  return a + (((float)rand())/RAND_MAX) * ( b - a );
}

void
createTest(int numLeafs,
	   int seqLen,
	   int ultrametricDeviation,
	   float diameterFactor,
	   char *seqfile_results,
	   char *treefile_results){
  
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
      int numWritten = sprintf(t,"%f",elen*
			       rand_float(((float)1)/ultrametricDeviation, 
					  ultrametricDeviation));
      t = t+numWritten;
    }
    else {
      *t = *s;
      t++; s++;
    }
  }
   //add the final semin colon (beep doesn't print it)
  *t = ';';
  t++;
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
  //execute seq-gen
  sprintf(str,  SEQGEN " -l%d "
          //" -wa " //write ancestor sequences  
          //JUKES CANTOR
          //	  " -s %f -n 1 -m HKY -t 0.5 -f  0.25 0.25 0.25 0.25 < "
          //" -n1 -mHKY -t0.5 -f0.25,0.25,0.25,0.25 < " //without diamter factor
          //K2P with fix ratio 2
          " -s%f -n1 -mHKY -t2 -f0.25,0.25,0.25,0.25 < "
          SEQGEN_IN " > " SEQGEN_OUT,
          seqLen, diameterFactor);
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

//    float diamterFactors[] = {0.05,0.1,0.25,0.5};
//    int numDF = 4;
//    int seqLengths[] = {50,100,250,500,1000,2000,4000};
//    int numSL = 6;
//    int ultrametricDeviation = 4;
//    int leafSizes[] = {10,25,50,100,200,400};
//    int numLS = 6;
//    int numTreesPerSize = 20;

  float diamterFactors[] = {0.05,0.25};
  int numDF = 2;
  int seqLengths[] = {50,100,250,500,1000,2000,4000};
  int numSL = 7;
  int ultrametricDeviation = 4;
  int leafSizes[] = {10,20,30,40,50,75};
  int numLS = 6;
  int numTreesPerSize = 20;

//   float diamterFactors[] = {0.25};
//   int numDF = 1;
//   int seqLengths[] = {500,1000,2000,4000};
//   int numSL = 4;
//   int ultrametricDeviation = 4;
//   int leafSizes[] = {10,20,30,40,50,60,70,80,90,100};
//   int numLS = 10;
//   int numTreesPerSize = 20;

//    float diamterFactors[] = {0.05};
//    int numDF = 1;
//    int seqLengths[] = {50};
//    int numSL = 1;
//    int ultrametricDeviation = 4;
//    int leafSizes[] = {10};
//    int numLS = 1;
//    int numTreesPerSize = 20;

  char treefilename[100];
  char seqfilename[100];
  //remove all test files
  system("rm -f " TESTDIR "test_leafs*" );

  int test = 1;
  //leaf sizes
  for ( int ls = 0 ; ls < numLS ; ls++ ){
    //diameter factor
    for ( int df = 0 ; df < numDF ; df++ ){
      //seq lengths
      for ( int sl = 0 ; sl < numSL ; sl++ ){
	//num trees for each combination
	sprintf(treefilename,TESTDIR "trees_leafs_%d_sqlen_%d_diamF_%f.txt",
		leafSizes[ls],seqLengths[sl],diamterFactors[df]);
	sprintf(seqfilename,TESTDIR "sequences_leafs_%d_sqlen_%d_diamF_%f.txt",
		leafSizes[ls],seqLengths[sl],diamterFactors[df]);
	//char numtests[100];
	//sprintf(numtests,"echo %d    num trees in file > %s", numTreesPerSize,filename);
	//system(numtests);
	for ( int tps = 0 ; tps < numTreesPerSize ; tps++ ){
	  cout <<"*************************** Doing test " << test << endl;
	  test++;
	
	  createTest(leafSizes[ls], seqLengths[sl], ultrametricDeviation,
		     diamterFactors[df], seqfilename,treefilename);
	}
      }
    }
  }

}




































