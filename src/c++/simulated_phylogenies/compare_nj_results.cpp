//--------------------------------------------------
//                                        
// File: compare_nj_results.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: compare_nj_results.cpp,v 1.2 2006/05/17 14:58:23 isaac Exp $                                 
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

typedef struct {
  int numLeafs;
  float diameterFactor;
  int seqLen;

  //Average RF
  float nj;
  float weighbor;
  float fnjcol;

} result;

vector<result> results;

float
getAverageRF(char *rfOutFile){

  ifstream tmp;
  open_read_stream(rfOutFile, tmp);
  char line1[100000];
  char line2[100000];
  char *line1ptr = line1;
  char *line2ptr = line2;
  char *tmpptr;  
  while ( !tmp.eof() ){
    tmpptr = line2ptr;
    line2ptr = line1ptr;
    line1ptr = tmpptr;
    tmp.getline(line1ptr,100000);
  }
  tmp.close();

  cout << line2ptr << endl;

  float avg;
  sscanf(line2ptr,"The total average RF dist was: %f",&avg);
  
  return avg;
}


void
runTest(int numLeafs,
        int seqLen,
        int ultrametricDeviation,
        float diameterFactor,
        char *seqfile_results,
        char *treefile_results){

  cout << "**************************" << endl;
  PRINT(numLeafs);PRINT(seqLen);PRINT(ultrametricDeviation);PRINT(diameterFactor);
  cout <<"---" << endl;

  //RUN NJ
  string njresult = string(seqfile_results) +".dist.njtree";
  string nj = "NJ " + string(seqfile_results) +".dist  > " + njresult;
  cout << nj << endl;
  system(nj.c_str());
  string rf_nj = "RF_dist " + string(treefile_results) + " " + njresult +" -m 20 > tmprf.txt";
  system(rf_nj.c_str());
  float nj_averageRF = getAverageRF("tmprf.txt");

  PRINT(nj_averageRF);

  //RUN NJ
  string wresult = string(seqfile_results) +".dist.wtree";
  
  string wstr = string("~/10giga_volume/Weighbor/weighbor -L ") +seqLen + " -i " + string(seqfile_results) +".dist  > " + wresult;
  cout << wstr << endl;
  system(wstr.c_str());
  string rf_w = "RF_dist " + string(treefile_results) + " " + wresult +" -m 20 > tmprf.txt";
  system(rf_w.c_str());
  float w_averageRF = getAverageRF("tmprf.txt");

  PRINT(w_averageRF);
  
  //RUN FNJ_collapse_method
  string fnjcollapse_result = string(seqfile_results) +".dist.fnjcollapse";
  string fnjcol = "FNJ_collapse_method -i " + string(seqfile_results) + ".dist -seqlen "+seqLen +" -m 20 > " + fnjcollapse_result;
  cout << fnjcol << endl;
  system(fnjcol.c_str());
  string rf_fnjcol = "RF_dist " + string(treefile_results) + " " + fnjcollapse_result +" -m 20 > tmprf.txt";
  system(rf_fnjcol.c_str());
  float fnjcol_averageRF = getAverageRF("tmprf.txt");

  PRINT(fnjcol_averageRF);


  result r;
  r.numLeafs = numLeafs;
  r.diameterFactor = diameterFactor;
  r.seqLen = seqLen;
  r.nj = nj_averageRF;
  r.weighbor = w_averageRF;
  r.fnjcol = fnjcol_averageRF;

  results.push_back(r);
}





int
main(int argc, char **argv){


  float diamterFactors[] = {0.05,0.1,0.25};//,0.5};
  int numDF = 3;
  int seqLengths[] = {50,100,250,500,1000,2000,4000};
  int numSL = 7;
  int ultrametricDeviation = 4;
  int leafSizes[] = {10,20,30,40,50,75,100,200,400};
  int numLS = 9;
  int numTreesPerSize = 20;

  int tot = numDF*numSL*numLS+2;

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
        runTest(leafSizes[ls], seqLengths[sl], ultrametricDeviation,
                diamterFactors[df], seqfilename,treefilename);

        if ( tot == 0 )
          goto end;
        tot--;
        cout << "******************* Tests left "<< tot << endl;
      }
    }
  }
 end:


  //MAKE PLOTS
  
  char datfilename[100];
  char plotfilename[100];
  for ( int df = 0 ; df < numDF ; df++ ){


    //---------------------------------------------------
    //MAKE PLOTS FOR DIFFERENT LEAFS SIZES
    for ( int ls = 0 ; ls < numLS ; ls++ ){
      
      sprintf(datfilename,"diamF_%f_leafs_%d.dat",
              diamterFactors[df],leafSizes[ls]);
      
      ofstream out;
      open_write_stream(datfilename, out);
      out << "# PLOT for fixed tree size and variable sequence length" << endl;
      out << "# Diamter Factor " << diamterFactors[df] << "     num leafs " << leafSizes[ls]  << endl;
      out <<"SeqLen \tNJ \tWeighbor \tFNJ_COL" << endl;

      for ( size_t i = 0 ; i < results.size() ; i++ ){
        result r = results[i];
        if ( r.diameterFactor != diamterFactors[df] )continue;
        if ( r.numLeafs != leafSizes[ls] )continue;
        out << r.seqLen << " \t" << r.nj << " \t" << r.weighbor  << " \t" << r.fnjcol <<endl;
      }
      out.close();
    
      sprintf(plotfilename,"diamF_%f_leafs_%d.eps",
              diamterFactors[df],leafSizes[ls]);
      cout << "***** DOing GNUPLOT" << endl;
      system("rm -f ddd.txt");
      open_write_stream("ddd.txt",out);
      out.precision(2);
      out << "reset\n"
          <<"set terminal postscript eps solid color \"Time-Roman\" 18\n"
          <<"set output '" << plotfilename <<"'\n"
          <<"set xlabel 'SeqLen'\n"
          <<"set ylabel 'Avg RF. " <<"'\n"
          <<"set title 'NJ,WeighBor, FNJ_COL  diamF "<< diamterFactors[df] <<", numLeafs. " << leafSizes[ls] << "'\n"
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
          <<"plot '"<< datfilename <<"' using 1:2 title 'NJ', '"<< datfilename <<"' using 1:3 title 'Weighbor', '" << datfilename <<"' using 1:4 title 'FNJ_COL'\n";

      out.close();
      string gnu("gnuplot ddd.txt");
      system(gnu.c_str());
    }
    //--------------------------------------------------
    //MAKE PLOTS FOR DIFFERENT SEQ LENGTHS
    for ( int sl = 0 ; sl < numSL ; sl++ ){
      sprintf(datfilename,"diamF_%f_seqLen_%d.dat",
              diamterFactors[df],seqLengths[sl]);
      
      ofstream out;
      open_write_stream(datfilename, out);
      out << "# PLOT for fixed sequence length and variable tree size" << endl;
      out << "# Diamter Factor " << diamterFactors[df] << "     sequence length " << seqLengths[sl]  << endl;
      out <<"NumLeafs \tNJ \tWeighbor \tFNJ_COL" << endl;

      for ( size_t i = 0 ; i < results.size() ; i++ ){
        result r = results[i];
        if ( r.diameterFactor != diamterFactors[df] )continue;
        if ( r.seqLen != seqLengths[sl] )continue;
        out << r.numLeafs << " \t" << r.nj << " \t" << r.weighbor << " \t" << r.fnjcol <<endl;
      }
      out.close();
      sprintf(plotfilename,"diamF_%f_seqLen_%d.eps",
              diamterFactors[df],seqLengths[sl]);
      cout << "***** DOing GNUPLOT" << endl;
      system("rm -f ddd.txt");
      open_write_stream("ddd.txt",out);
      out.precision(2);
      out << "reset\n"
          <<"set terminal postscript eps solid color \"Time-Roman\" 18\n"
          <<"set output '" << plotfilename <<"'\n"
          <<"set xlabel 'Num Leafs'\n"
          <<"set ylabel 'Avg RF. " <<"'\n"
          <<"set title 'NJ,Weighbor,FNJ_COLLAPSE  diamF "<< diamterFactors[df] <<", seqLen. " << seqLengths[sl] << "'\n"
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
          <<"plot '"<< datfilename <<"' using 1:2 title 'NJ', '"<< datfilename <<"' using 1:3 title 'Weighbor', '" << datfilename <<"' using 1:4 title 'FNJ_COL'\n";

      out.close();
      string gnu("gnuplot ddd.txt");
      system(gnu.c_str());
    }
  }
  return 1;
}
