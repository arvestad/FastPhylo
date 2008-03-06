//--------------------------------------------------
//                                        
// File: BootstrapTest.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: BootstrapTest.cpp,v 1.12 2006/09/04 11:50:12 isaac Exp $
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
#include "SequenceTree.hpp"
#include "Sequence.hpp"
#include "SequenceParsimony.hpp"
#include "LeastSquaresFit.hpp"

#include <vector>
#include "arg_utils.h"

#define TESTDIR "testfiles/"

bool rerun_algos = true;

using namespace std;

typedef struct {
  int numLeafs;
  float diameterFactor;
  int seqLen;

  //L_2 fit
  float avg_model_l2fit;
  float avg_rf_best_l2fit;
  float avg_l2fit;

  //Average Parsimony difference compared to model tree
  float nj_avg_parsi_diff;
  float best_parsi_avg_diff;
  float avg_parsi_diff;
  //Average RF
  float best_parsi;
  float nj;
  float best_boot;
  float worst_boot;
  float average_boot;
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
runTest(int num_datasets,
	int num_boot_each_dataset,
	int numLeafs,
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

  //The file is divided into 20 different datasets each of which have been bootstrapped
  //1000 times in addition to the original sequence. Thus there are 1001*20 matrices.  
  
  //RUN NJ
  //compute all nj trees 
  string njresult = string(seqfile_results) +".dist.njtree";
  string nj = "NJ " + string(seqfile_results) +".dist  > " + njresult;
  cout << nj << endl;
  t = time(&t);
  if(rerun_algos){
    system(nj.c_str());
    cout << "TIME: " << difftime(time(NULL),t) << "sec" << endl;
  }

  ifstream modelTreeFile;
  open_read_stream(treefile_results,modelTreeFile);
  ifstream njTreeFile;
  open_read_stream(njresult.c_str(),njTreeFile);
  ifstream seqFile;
  open_read_stream(seqfile_results, seqFile);
  ifstream dmFile;
  open_read_stream((string(seqfile_results) +".dist").c_str(),dmFile);
  StrDblMatrix tmpdm(100);
  StrDblMatrix origdm(100);

  double sum_model_l2fit=0;
  double sum_avg_all_l2fit=0;
  double sum_best_l2fit=0;
  double sum_rf_ofbestl2=0;

  double sum_nj_parsi_diff=0;
  double sum_best_parsi_diff=0;
  double sum_avg_parsi_diff=0;

  double sum_best_parsi=0;
  double sum_nj=0;
  double sum_best_boot=0;
  double sum_worst_boot=0;
  double sum_average_boot=0;
  for(int dataset=0 ; dataset<num_datasets ; dataset++ ){
    SEPARATOR();PRINT(dataset);
    //read model tree    
    SequenceTree modelTree(modelTreeFile);
    PRINT(modelTree);

    vector<Sequence> seqs;
    Sequence::readSequences(seqs, seqFile);
    mapSequencesOntoTree(modelTree,seqs);
    shortcutDegree2Nodes(modelTree);
    size_t modelP = computeMostParsimoniousSequences(modelTree);
    PRINT(modelP);

    //read the original nj tree
    SequenceTree tmpTree(njTreeFile);
    double rf_orig=computeRobinsonFoulds(tmpTree,modelTree);
    PRINT(rf_orig);
    mapSequencesOntoTree(tmpTree,seqs);
    size_t njP = computeMostParsimoniousSequences(tmpTree);
    PRINT(njP);

    //read the original distance matrix
    origdm.fillFromStream(dmFile); 
    double modelFit = computeLeastSquaresEdgeLengths(origdm,modelTree);

    double njFit = computeLeastSquaresEdgeLengths(origdm,tmpTree);
    double all_l2fit = njFit;
    double bestl2fit = njFit;
    double rf_ofbestl2 = rf_orig;
    PRINT(modelFit);PRINT(njFit);

    //read the 1000 boot trees
    size_t best_parsi_score=njP;
    double sum_parsi_diff=(njP-modelP);
    double rf_best_parsi = rf_orig;
    double best_boot=rf_orig;
    double worst_boot=rf_orig;
    double sum_boot=rf_orig;
    for(int boot=0 ; boot<num_boot_each_dataset; boot++){
      tmpdm.fillFromStream(dmFile);
      tmpTree = SequenceTree(njTreeFile);
      double rf_boot=computeRobinsonFoulds(tmpTree,modelTree);
      mapSequencesOntoTree(tmpTree,seqs);
      size_t bootP = computeMostParsimoniousSequences(tmpTree);
      
      double bootFit = computeLeastSquaresEdgeLengths(origdm,tmpTree);
      all_l2fit+=bootFit;
      PRINT(bootFit);
      if(bootFit<bestl2fit){
	bestl2fit=bootFit;
	rf_ofbestl2 = rf_boot;
      }

      PRINT(bootP);
      sum_parsi_diff += (bootP-modelP);
      if(bootP<best_parsi_score){
	rf_best_parsi = rf_boot;
	best_parsi_score = bootP;
      }
      PRINT(boot);PRINT(rf_boot);
      if(rf_boot<best_boot)
	best_boot = rf_boot;
      if(rf_boot>worst_boot)
	worst_boot = rf_boot;
      sum_boot+=rf_boot;
    }

    SEPARATOR();PRINT(rf_orig);PRINT(best_boot);PRINT(worst_boot);PRINT(sum_boot/(1+num_boot_each_dataset));
    PRINT(modelP);PRINT(njP);PRINT(best_parsi_score);PRINT(rf_best_parsi);PRINT(modelFit);PRINT(bestl2fit);PRINT(rf_ofbestl2);

    sum_model_l2fit += modelFit;
    sum_best_l2fit += bestl2fit;
    sum_rf_ofbestl2 += rf_ofbestl2;
    sum_avg_all_l2fit += all_l2fit/(1+num_boot_each_dataset);

    sum_nj_parsi_diff += (njP-modelP);
    sum_best_parsi_diff += (best_parsi_score-modelP);
    sum_avg_parsi_diff += sum_parsi_diff/(1+num_boot_each_dataset);
    sum_best_parsi += rf_best_parsi;
    sum_nj += rf_orig;
    sum_best_boot += best_boot;
    sum_worst_boot +=worst_boot;
    sum_average_boot += sum_boot/(1+num_boot_each_dataset);//including the original nj
  }
  
  //the result is the average over 20 runs
  r.avg_l2fit = sum_avg_all_l2fit/num_datasets;
  r.avg_rf_best_l2fit = sum_rf_ofbestl2/num_datasets;
  r.avg_model_l2fit = sum_model_l2fit/num_datasets;

  r.nj_avg_parsi_diff = sum_nj_parsi_diff/num_datasets;
  r.best_parsi_avg_diff = sum_best_parsi_diff/num_datasets;
  r.avg_parsi_diff = sum_avg_parsi_diff/num_datasets;
  r.best_parsi = sum_best_parsi/num_datasets;
  r.nj = sum_nj/num_datasets;
  r.best_boot = sum_best_boot/num_datasets;
  r.worst_boot = sum_worst_boot/num_datasets;
  r.average_boot = sum_average_boot/num_datasets;
 

  //****
  modelTreeFile.close();
  njTreeFile.close();
  seqFile.close();
  dmFile.close();
  //****
  results.push_back(r);
}





int
main(int argc, char **argv){

  if(HAS_OPTION("-no_rerun"))rerun_algos = false;

  float diamterFactors[] = {0.05,0.25};
  int numDF = 2;
  int seqLengths[] = {50,100,250,500,1000,2000,4000};
  int numSL = 5;
  int ultrametricDeviation = 4;
  int leafSizes[] = {10,20,30,40,50,75};
  int numLS = 3;
  int numTreesPerSize = 20;
  int numBootsPerTree = 100;

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
        runTest(numTreesPerSize, numBootsPerTree, leafSizes[ls], seqLengths[sl], ultrametricDeviation,
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
  latex <<"\\documentclass{article}\\usepackage{epic}\\usepackage{graphics}"<< endl
	<<"\\textwidth 6.8in \\textheight 9.275in \\topmargin -0.475in \\oddsidemargin -.25in \\evensidemargin -.25in"<< endl
	<< "\\begin{document}"<<endl;
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
            << "\\end{minipage}\\begin{minipage}{0.9\\textwidth}" << endl
            << "\\begin{tabular}{r|rrrrrr}\\textbf{Length} & \\textbf{NJ} & \\textbf{Best} & \\textbf{Avg.} & \\textbf{Worst}& \\textbf{BestP}&\\textbf{BestL2}\\\\" << endl << "\\hline" << endl;
      
      ofstream out;
      open_write_stream(datfilename, out);
      out << "# PLOT for fixed tree size and variable sequence length" << endl
          << "# Diamter Factor " << diamterFactors[df] << "     num leafs " << leafSizes[ls]  << endl
          <<"SeqLen \tNJ \tBest \tAvg. \tWorst \tBestParsi \tBestL2" << endl;
      for ( size_t i = 0 ; i < results.size() ; i++ ){
        result r = results[i];
        if ( r.diameterFactor != diamterFactors[df] )continue;
        if ( r.numLeafs != leafSizes[ls] )continue;
	latex << r.seqLen << " & " << r.nj <<"  & " << r.best_boot << " & " << r.average_boot  
	      << " & " << r.worst_boot << " & " << r.best_parsi << " & " << r.avg_rf_best_l2fit<<"\\\\" <<endl;
	out << r.seqLen << " \t" << r.nj <<" \t" << r.best_boot << " \t" << r.average_boot  
	    << " \t" << r.worst_boot << " \t " << r. best_parsi << " \t " << r.avg_rf_best_l2fit << endl;
        
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

      out <<"plot '"<< datfilename <<"' using 1:2 title 'NJ', '"<< datfilename <<"' using 1:3 title 'Best', '"
          << datfilename <<"' using 1:4 title 'Avg.', '" << datfilename <<
	"' using 1:5 title 'Worst', '" << datfilename <<"' using 1:6 title 'BestP', '" << datfilename <<"' using 1:7 title 'BestL2'\n";
	
      out.close();
      string gnu("gnuplot ddd.txt");
      system(gnu.c_str());
    }
    //latex << ""<<endl;

    latex << "\\caption{The number of leafs is fixed to " << leafSizes[ls] << " and the length of the sequences varies for different diameter factors.}\\end{figure}" << endl<<endl;
    latex <<"$\\;$\\newline"<<endl;
  }

  //-------------------------------
  latex << endl<< endl<< endl<< endl<<
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< endl<<
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< endl<<
    "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<< endl;
  latex << "\\section{Fixed Sequence Length - Varaying Tree Size}" << endl << endl<<endl;

  
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
            << "\\end{minipage}\\begin{minipage}{0.9\\textwidth}" << endl
            << "\\begin{tabular}{r|rrrrr}\\textbf{Leafs} & \\textbf{NJ} & \\textbf{Best} & \\textbf{Avg.} & \\textbf{Worst} & \\textbf{BestP}\\\\" << endl << "\\hline" << endl;
      
      ofstream out;
      open_write_stream(datfilename, out);
      out << "# PLOT for fixed sequence length and variable tree size" << endl
          << "# Diamter Factor " << diamterFactors[df] << "     sequence length " << seqLengths[sl]  << endl
          <<"Leafs \tNJ \tBest \tAvg. \tWorst \tBestP" << endl;
      for ( size_t i = 0 ; i < results.size() ; i++ ){
        result r = results[i];
        if ( r.diameterFactor != diamterFactors[df] )continue;
        if ( r.seqLen != seqLengths[sl] )continue;
        
	latex << r.numLeafs << " & " << r.nj <<"  & " << r.best_boot << " & " << r.average_boot  
	      << " & " << r.worst_boot<< " & " << r.best_parsi << "\\\\" <<endl;
	out << r.numLeafs << " \t" << r.nj <<" \t" << r.best_boot << " \t" << r.average_boot  
	    << " \t" << r.worst_boot << " \t" << r.best_parsi <<endl;
        
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
          <<"plot '"<< datfilename <<"' using 1:2 title 'NJ', '"<< datfilename <<"' using 1:3 title 'Best', '"
          << datfilename <<"' using 1:4 title 'Avg.', '" << datfilename 
	  <<"' using 1:5 title 'Worst', '" << datfilename <<"' using 1:6 title 'BestP'\n";

      out.close();
      string gnu("gnuplot ddd.txt");
      system(gnu.c_str());
    }
    latex <<"\\caption{The sequence length is fixed to " << seqLengths[sl] << " and the number of trees varies for different diameter factors.}\\end{figure}" << endl<<endl;
    latex <<"$\\;$\\newline"<<endl;
  }

  //---------------------------------------------
  latex << endl<< endl<< endl<< endl;
  latex << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  latex << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  latex << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
  latex << "\\section{Fixed Tree Size - Varaying Sequence Length PARSIMONY}" << endl << endl<<endl;

  for ( int ls = 0 ; ls < numLS ; ls++ ){
    latex << "%%%%%%%" << endl;
    //latex << "\\subsection{Trees with " << leafSizes[ls]<<" leafs}"<< endl;
    latex << "\\begin{figure}[h]\\centering"<<endl;
    for ( int df = 0 ; df < numDF ; df++ ){
      sprintf(datfilename,"diamF_%f_seqlen_%dPARSI.dat",diamterFactors[df],leafSizes[ls]);
      sprintf(plotfilename,"diamF_%f_seqlen_%dPARSI.eps",diamterFactors[df],leafSizes[ls]);

      //--
      latex << "%Diamter Factor " << diamterFactors[df] << "     num leafs " << leafSizes[ls]  << endl
	    << "% THE PARSIMONY DIFFERENCE " << endl
            << "\\begin{minipage}{0.45\\textwidth}"<< endl
            << "\\resizebox{55mm}{!}{\\includegraphics{" << plotfilename << "}}" <<endl
            << "\\end{minipage}\\begin{minipage}{0.9\\textwidth}" << endl
            << "\\begin{tabular}{r|rrr}\\textbf{Length} & \\textbf{NJ} & \\textbf{Best} & \\textbf{Avg.} \\\\" << endl << "\\hline" << endl;
      
      ofstream out;
      open_write_stream(datfilename, out);
      out << "# PLOT for fixed tree size and variable sequence length" << endl
          << "# Diamter Factor " << diamterFactors[df] << "     num leafs " << leafSizes[ls]  << endl
	  << "# PARSIMONY DIFFERENCE" <<endl
          <<"SeqLen \tNJ \tBest \tAvg." << endl;
      for ( size_t i = 0 ; i < results.size() ; i++ ){
        result r = results[i];
        if ( r.diameterFactor != diamterFactors[df] )continue;
        if ( r.numLeafs != leafSizes[ls] )continue;
	latex << r.seqLen << " & " << r.nj_avg_parsi_diff <<"  & " << r.best_parsi_avg_diff << " & " << r.avg_parsi_diff
	      << "\\\\" <<endl;
	out << r.seqLen << " \t" << r.nj_avg_parsi_diff <<" \t" << r.best_parsi_avg_diff << " \t" << r.avg_parsi_diff << endl;
        
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
          <<"set ylabel 'Avg Parsi Diff. " <<"'\n"
          <<"set title 'Diameter factor = "<< diamterFactors[df] <<", Number of leafs = " << leafSizes[ls] << "'\n"
          <<"set auto y\n"
          <<"set style data line\n"
          <<"set key left bottom\n"
          <<"set style fill solid border -1\n"
          <<"#set xtic rotate by 90\n"
          <<"#set bmargin 10\n";

      out <<"plot '"<< datfilename <<"' using 1:2 title 'NJ-diff', '"<< datfilename <<"' using 1:3 title 'Best diff', '"
          << datfilename <<"' using 1:4 title 'Avg. diff'\n";
	
      out.close();
      string gnu("gnuplot ddd.txt");
      system(gnu.c_str());
    }
    //latex << ""<<endl;

    latex << "\\caption{The number of leafs is fixed to " << leafSizes[ls] << " and the length of the sequences varies for different diameter factors.}\\end{figure}" << endl<<endl;
    latex <<"$\\;$\\newline"<<endl;
  }

  
  latex << endl<< endl<< endl<< endl<<
  
  //-------------------------------------------------
  latex <<"\\end{document}"<< endl;
  latex.close();
  
  return 1;
}

