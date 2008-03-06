//--------------------------------------------------
//                                        
// File: compare_aml_results.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: compare_aml_results.cpp,v 1.2 2006/03/27 14:32:27 isaac Exp $                                 
//
//--------------------------------------------------
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

#include <vector>

#define TESTDIR "testfiles/"

using namespace std;

void
readLiklihoodFloats(vector<double> &dvec){

  ifstream tmp;
  open_read_stream("tmp.txt", tmp);
  char line[1000];
  
  while ( !tmp.eof() ){
    tmp.getline(line,1000);
    float f;
    if ( sscanf(line, "Loglikelihood = %f",&f) != 1 )
      continue;

    dvec.push_back((double)f);
  }

  tmp.close();

}

double
normalizedDiff(vector<double> &a,
               vector<double> &b){
  assert ( a.size() == b.size() );
  double tmp = 0;
  for ( size_t i = 0 ; i < a.size() ; i++ ){
    tmp += ( a[i]-b[i] )/b[i];
  }

  tmp = tmp / a.size();
  return tmp;
}

void
compareAlgoTest(char *seqfile,
                char *treefile,
                double &bigdiff,
                double &liftdiff,
                double &parsdiff){

  string amlbasestr("aml_tree -seqgen -d 20 -seqs ");
  amlbasestr += seqfile;

  //execute big AML
  string big = amlbasestr + "  > tmp.txt";
  system(big.c_str());
  vector<double> bigvec;
  readLiklihoodFloats(bigvec);

  //execute leaf lift AML
  string lift = amlbasestr + " -tree " + treefile + "  > tmp.txt";
  system(lift.c_str());
  vector<double> liftvec;
  readLiklihoodFloats(liftvec);

  //execute parsimony AML
  lift = amlbasestr + " -tree " + treefile + " -use-parsimony  > tmp.txt";
  system(lift.c_str());
  vector<double> parsvec;
  readLiklihoodFloats(parsvec);

  //execute ancestoral AML
  lift = amlbasestr + " -tree " + treefile + " -use-ancestral  > tmp.txt";
  system(lift.c_str());
  vector<double> avec;
  readLiklihoodFloats(avec);


  bigdiff = normalizedDiff(bigvec,avec);
  liftdiff = normalizedDiff(liftvec,avec);
  parsdiff = normalizedDiff(parsvec,avec);
}


int
main(int argc, char **argv){

  cout << "run as: ./compare_aml_results 0.25\nThe argument is the mutation probability for branches" << endl;
  float mutationProb = atof(argv[1]);
  int seqLengths[] = {10,50,100,500,1000,6000};
  int numSL = 6;
  int leafSizes[] = {5,10,15,20,30,40,50,60,70,80,90,100};
  int numLS = 12;


  char treefilename[100];
  char seqfilename[100];



  for ( int leni = 0 ; leni < numSL ; leni++ ){

    vector<double> big(numLS);
    vector<double> lift(numLS);
    vector<double> pars(numLS);
    for ( int li = 0 ; li < numLS ; li++ ){
      sprintf(treefilename,TESTDIR "trees_leafs_%d_sqlen_%d_mprob_%f.txt",
              leafSizes[li],seqLengths[leni],mutationProb);
      sprintf(seqfilename,TESTDIR "sequences_leafs_%d_sqlen_%d_mprob_%f.txt",
              leafSizes[li],seqLengths[leni],mutationProb);

      PRINT(seqfilename);
      double bigdiff;
      double liftdiff;
      double pdiff;
      compareAlgoTest(seqfilename,treefilename,bigdiff,liftdiff,pdiff);
      big[li] = bigdiff;
      lift[li] = liftdiff;
      pars[li] = pdiff;
      PRINT(bigdiff);PRINT(liftdiff);PRINT(pdiff);
    }

    ofstream out;
    string fn("comparealgos_seqlen_");
    fn = fn + seqLengths[leni];
    fn += "mutationProb_0_";
    fn = fn+mutationProb+".dat";
    open_write_stream(fn.c_str(),out);
    out <<"****************************"<<endl;
    out <<"#The normalized difference between the original\n#ancestoral sequences and the output of the resp algos."<< endl;
    out <<"#Each algorithm (and the orignal ancestoral) sequences have been improved by the local improvement heursitic" << endl;
    out <<"#Sequence length " << seqLengths[leni] << endl;
    out <<"Leafs    Big_AML   Small_AML    Parsimony" << endl; 
    for ( int li = 0 ; li < numLS ; li++ ){
      out << leafSizes[li] << "   " << big[li] << "   " << lift[li] << "   " << pars[li] << endl;
    }
    out.close();
    string cat = "cat " + fn;
    system(cat.c_str());

    cout << "***** DOing GNUPLOT" << endl;
    string gn("plot_seqlen_");
    gn = gn+ seqLengths[leni];
    gn = gn + "_mProb_" + mutationProb;
    gn += ".eps";
    system("rm -f ddd.txt");
    open_write_stream("ddd.txt",out);
    out.precision(2);
    out << "reset\n"
        <<"set terminal postscript eps solid color \"Time-Roman\" 18\n"
        <<"set output '" << gn <<"'\n"
        <<"set xlabel 'Num leafs'\n"
        <<"set ylabel 'Norm diff. " <<"'\n"
        <<"set title 'AML algorithms, m.prob. 0-"<< mutationProb <<", seq.len. " << seqLengths[leni] << "'\n"
      //<<"set datafile missing \"-\"\n"
      //<<"set xrange [-1:10]\n"
        <<"set auto y\n"
      //       <<"set style data histogram\n"
        <<"set style data line\n"
        <<"set key left bottom\n"
      //<<"set style histogram cluster gap 1\n"
        <<"set style fill solid border -1\n"
      //        <<"set xtics border (\"10\" 0.00000, \"20\" 1.00000, \"30\" 2.00000, \"40\" 3.00000, \"50\" 4.00000, \"60\" 5.00000, \"70\" 6.00000, \"80\" 7.00000, \"90\" 8.00000, \"100\" 9.00000)\n"
      //        <<"set boxwidth 0.9\n"
        <<"#set xtic rotate by 90\n"
        <<"#set bmargin 10\n"
        <<"plot '"<< fn <<"' using 1:2 title 'Big AML', '" << fn <<"' using 1:3 title 'Small AML',  '" << fn <<"' using 1:4 title 'Parsimony'\n";

      out.close();
      string gnu("gnuplot ddd.txt");
      system(gnu.c_str());
  }
}

