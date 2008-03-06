//--------------------------------------------------
//                                        
// File: compare_nj_results.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: compare_nj_bionj_weighbor_fnjcol.cpp,v 1.2 2006/05/22 13:58:23 isaac Exp $                                 
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
#include <time.h>

#include <vector>
#include "arg_utils.h"

#define TESTDIR "testfiles/"

bool rerun_algos = true;

using namespace std;

typedef struct {
  int numLeafs;
  float diameterFactor;
  int seqLen;

  //Average RF
  float nj;
  float bionj;
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

  result r;
  r.numLeafs = numLeafs;
  r.diameterFactor = diameterFactor;
  r.seqLen = seqLen;
  
  time_t t;
  
  
  //RUN NJ
  string njresult = string(seqfile_results) +".dist.njtree";
  string nj = "NJ " + string(seqfile_results) +".dist  > " + njresult;
  cout << nj << endl;
  t = time(&t);
  if(rerun_algos){
    system(nj.c_str());
    cout << "TIME: " << difftime(time(NULL),t) << "sec" << endl;
  }
  string rf_nj = "RF_dist " + string(treefile_results) + " " + njresult +" -m 20 > tmprf.txt";
  system(rf_nj.c_str());
  float nj_averageRF = getAverageRF("tmprf.txt");

  PRINT(nj_averageRF);
  r.nj = nj_averageRF;

  //RUN BIONJ
  string bionjresult = string(seqfile_results) +".dist.bionjtree";
  string bionj = "BIONJ " + string(seqfile_results) +".dist  > " + bionjresult;
  cout << bionj << endl;
  t = time(NULL);
  if(rerun_algos){
    system(bionj.c_str());
    cout << "TIME: " <<  difftime(time(NULL),t) << "sec" << endl;
  }
  string rf_bionj = "RF_dist " + string(treefile_results) + " " + bionjresult +" -m 20 > tmprf.txt";
  system(rf_bionj.c_str());
  float bionj_averageRF = getAverageRF("tmprf.txt");

  PRINT(bionj_averageRF);
  r.bionj = bionj_averageRF;
  //RUN WEIGHBOR
  if ( numLeafs <= 200 ){
    string wresult = string(seqfile_results) +".dist.wtree";

    string wstr = string("~/10giga_volume/Weighbor/weighbor -L ") +seqLen + " -i " + string(seqfile_results) +".dist  > " + wresult;
    cout << wstr << endl;
    t = time(NULL);
    if(rerun_algos){
      system(wstr.c_str());
      cout << "TIME: " <<  difftime(time(NULL),t) << "sec" << endl;
    }
    string rf_w = "RF_dist " + string(treefile_results) + " " + wresult +" -m 20 > tmprf.txt";
    system(rf_w.c_str());
    float w_averageRF = getAverageRF("tmprf.txt");

    PRINT(w_averageRF);
    r.weighbor = w_averageRF;
  }
  else{
    r.weighbor = -1;
    cout << "***More than 400 leafs... not running weigbor" <<endl;
  }

  //RUN FNJ_collapse_method
  string fnjcollapse_result = string(seqfile_results) +".dist.fnjcollapse";
  string fnjcol = "FNJ_collapse_method -i " + string(seqfile_results) + ".dist -seqlen "+seqLen +" -m 20 > " + fnjcollapse_result;
  cout << fnjcol << endl;
  t = time(NULL);
  if(rerun_algos){
    system(fnjcol.c_str());
    cout << "TIME: " << difftime(time(NULL),t) << "sec" << endl;
  }
  string rf_fnjcol = "RF_dist " + string(treefile_results) + " " + fnjcollapse_result +" -m 20 > tmprf.txt";
  system(rf_fnjcol.c_str());
  float fnjcol_averageRF = getAverageRF("tmprf.txt");

  PRINT(fnjcol_averageRF);
  r.fnjcol = fnjcol_averageRF;


  //****
  results.push_back(r);
}





int
main(int argc, char **argv){

  if(HAS_OPTION("-no_rerun"))rerun_algos = false;

  float diamterFactors[] = {0.05,0.1,0.25,0.5};
  int numDF = 4;
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


  //---------------------------------------------------------

  //MAKE PLOTS
  
  char datfilename[100];
  char plotfilename[100];
  ofstream latex;
  open_write_stream("Tests_on_simulated_data.tex", latex);
  latex <<"\\documentclass{article}\\usepackage{epic}\\usepackage{graphics}\\begin{document}"<<endl;
  latex << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  latex << "\\section{Fixed Tree Size - Varaying Sequence Length}" << endl << endl<<endl;

  for ( int ls = 0 ; ls < numLS ; ls++ ){
    latex << "%%%%%%%" << endl;
    //latex << "\\subsection{Trees with " << leafSizes[ls]<<" leafs}"<< endl;
    latex << "\\begin{figure}[h]\\centering"<<endl;
    for ( int df = 0 ; df < numDF ; df++ ){
      sprintf(datfilename,"diamF_%f_leafs_%d.dat",diamterFactors[df],leafSizes[ls]);
      sprintf(plotfilename,"diamF_%f_leafs_%d.eps",diamterFactors[df],leafSizes[ls]);

      //--
      latex << "%Diamter Factor " << diamterFactors[df] << "     num leafs " << leafSizes[ls]  << endl
            << "\\begin{minipage}{0.45\\textwidth}"<< endl
            << "\\resizebox{55mm}{!}{\\includegraphics{" << plotfilename << "}}" <<endl
            << "\\end{minipage}\\begin{minipage}{0.45\\textwidth}" << endl
            << "\\begin{tabular}{r|rrrr}\\textbf{Length} & \\textbf{NJ} & \\textbf{BioNJ} & \\textbf{Weigh.} & \\textbf{Coll.}\\\\" << endl << "\\hline" << endl;
      
      ofstream out;
      open_write_stream(datfilename, out);
      out << "# PLOT for fixed tree size and variable sequence length" << endl
          << "# Diamter Factor " << diamterFactors[df] << "     num leafs " << leafSizes[ls]  << endl
          <<"SeqLen \tNJ \tBioNJ \tWeighbor \tCOLL" << endl;
      for ( size_t i = 0 ; i < results.size() ; i++ ){
        result r = results[i];
        if ( r.diameterFactor != diamterFactors[df] )continue;
        if ( r.numLeafs != leafSizes[ls] )continue;
        if ( r.weighbor >= 0 ){
          latex << r.seqLen << " & " << r.nj <<"  & " << r.bionj << " & " << r.weighbor  << " & " << r.fnjcol << "\\\\" <<endl;
          out << r.seqLen << " \t" << r.nj <<" \t" << r.bionj << " \t" << r.weighbor  << " \t" << r.fnjcol <<endl;
        }
        else{
          latex << r.seqLen << " & " << r.nj <<"  & " << r.bionj << " & - & " << r.fnjcol << "\\\\" <<endl;
          out << r.seqLen << " \t" << r.nj <<" \t" << r.bionj << " \t - \t" << r.fnjcol <<endl;
        }
      }
      out.close();
      latex << "\\end{tabular}\\end{minipage}"<<endl;


      //make the eps
      cout << "***** DOing GNUPLOT" << endl;
      system("rm -f ddd.txt");
      open_write_stream("ddd.txt",out);
      out.precision(2);
      out << "reset\n"
          <<"set terminal postscript eps solid color \"Time-Roman\" 18\n"
          <<"set output '" << plotfilename <<"'\n"
          <<"set xlabel 'SeqLen'\n"
          <<"set ylabel 'Avg RF. " <<"'\n"
          <<"set title 'Diameter factor = "<< diamterFactors[df] <<", Number of leafs = " << leafSizes[ls] << "'\n"
          <<"set auto y\n"
          <<"set style data line\n"
          <<"set key left bottom\n"
          <<"set style fill solid border -1\n"
          <<"#set xtic rotate by 90\n"
          <<"#set bmargin 10\n";
      if(leafSizes[ls]<=200){
        out <<"plot '"<< datfilename <<"' using 1:2 title 'NJ', '"<< datfilename <<"' using 1:3 title 'BioNJ', '"
          << datfilename <<"' using 1:4 title 'Weighbor', '" << datfilename <<"' using 1:5 title 'Collapse'\n";
      }
      else{
        out <<"plot '"<< datfilename <<"' using 1:2 title 'NJ', '"<< datfilename <<"' using 1:3 title 'BioNJ', '"
            << datfilename <<"' using 1:5 title 'Collapse'\n";
      }
      out.close();
      string gnu("gnuplot ddd.txt");
      system(gnu.c_str());
    }
    //latex << ""<<endl;

    latex << "\\caption{The number of leafs is fixed to " << leafSizes[ls] << " and the length of the sequences varies for different diameter factors.}\\end{figure}" << endl<<endl;
    latex <<"$\\;$\\newline"<<endl;
  }

  latex << endl<< endl<< endl<< endl<<
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< endl<<
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< endl<<
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< endl;
  latex << "\\section{Fixed Sequence Length - Varaying Tree Size}" << endl << endl<<endl;
  //-----------------------
  for ( int sl = 0 ; sl < numSL ; sl++ ){
    latex << "%%%%%%%" << endl;
    latex << "\\begin{figure}[h]\\centering" <<endl;
    //latex << "\\subsection{Sequences of length " << seqLengths[sl]<<"}"<< endl;
    //latex << "The length of the sequences is fixed and the number leasf in the trees varies for different diameter factors."<<endl;

    for ( int df = 0 ; df < numDF ; df++ ){      
      sprintf(plotfilename,"diamF_%f_seqLen_%d.eps",diamterFactors[df],seqLengths[sl]);
      sprintf(datfilename,"diamF_%f_seqLen_%d.dat",diamterFactors[df],seqLengths[sl]);
    

      //--
      latex << "%Diamter Factor " << diamterFactors[df] << "   Sequence Length " << seqLengths[sl]  << endl
            << "\\begin{minipage}{0.45\\textwidth}"<< endl
            << "\\resizebox{55mm}{!}{\\includegraphics{" << plotfilename << "}}" <<endl
            << "\\end{minipage}\\begin{minipage}{0.45\\textwidth}" << endl
            << "\\begin{tabular}{r|rrrr}\\textbf{Leafs} & \\textbf{NJ} & \\textbf{BioNJ} & \\textbf{Weigh.} & \\textbf{Coll.}\\\\" << endl << "\\hline" << endl;
      
      ofstream out;
      open_write_stream(datfilename, out);
      out << "# PLOT for fixed sequence length and variable tree size" << endl
          << "# Diamter Factor " << diamterFactors[df] << "     sequence length " << seqLengths[sl]  << endl
          <<"Leafs \tNJ \tBioNJ \tWeighbor \tCOLL" << endl;
      for ( size_t i = 0 ; i < results.size() ; i++ ){
        result r = results[i];
        if ( r.diameterFactor != diamterFactors[df] )continue;
        if ( r.seqLen != seqLengths[sl] )continue;
        if ( r.weighbor >= 0 ){
          latex << r.numLeafs << " & " << r.nj <<"  & " << r.bionj << " & " << r.weighbor  << " & " << r.fnjcol << "\\\\" <<endl;
          out << r.numLeafs << " \t" << r.nj <<" \t" << r.bionj << " \t" << r.weighbor  << " \t" << r.fnjcol <<endl;
        }
        else{
          latex << r.numLeafs << " & " << r.nj <<"  & " << r.bionj << " & - & " << r.fnjcol << "\\\\" <<endl;
          out << r.numLeafs << " \t" << r.nj <<" \t" << r.bionj << " \t - \t" << r.fnjcol <<endl;
        }

      }
      out.close();
      latex << "\\end{tabular}\\end{minipage}"<<endl ;

      //make the eps
      cout << "***** DOing GNUPLOT" << endl;
      system("rm -f ddd.txt");
      open_write_stream("ddd.txt",out);
      out.precision(2);
      out << "reset\n"
          <<"set terminal postscript eps solid color \"Time-Roman\" 18\n"
          <<"set output '" << plotfilename <<"'\n"
          <<"set xlabel 'Num. Leafs.'\n"
          <<"set ylabel 'Avg RF. " <<"'\n"
          <<"set title 'Diameter factor = "<< diamterFactors[df] <<", Sequences of length = " << seqLengths[sl] << "'\n"
          <<"set auto y\n"
          <<"set style data line\n"
          <<"set key left bottom\n"
          <<"set style fill solid border -1\n"
          <<"#set xtic rotate by 90\n"
          <<"#set bmargin 10\n"
          <<"plot '"<< datfilename <<"' using 1:2 title 'NJ', '"<< datfilename <<"' using 1:3 title 'BioNJ', '"
          << datfilename <<"' using 1:4 title 'Weighbor', '" << datfilename <<"' using 1:5 title 'Collapse'\n";

      out.close();
      string gnu("gnuplot ddd.txt");
      system(gnu.c_str());
    }
    latex <<"\\caption{The sequence length is fixed to " << seqLengths[sl] << " and the number of trees varies for different diameter factors.}\\end{figure}" << endl<<endl;
    latex <<"$\\;$\\newline"<<endl;
  }

  latex <<"\\end{document}"<< endl;
  latex.close();
  
  return 1;
}
