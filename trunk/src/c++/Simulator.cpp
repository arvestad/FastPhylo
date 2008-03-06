
#include "Simulator.hpp"
#include "std_c_utils.h"
#include "log_utils.hpp"
#include <sstream>
#include <math.h>

#define BEEP "~/10giga_volume/beep/beep_generateTree"
#define BEEP_OUT "beep_outfile.txt"
#define SEQGEN "~/10giga_volume/Seq-Gen.v1.3.2/source/seq-gen"
#define SEQGEN_IN "seqgen_infile.txt"
#define SEQGEN_OUT "seqgen_outfile.txt"


#define ROSE "~/10giga_volume/rose-1.3/src/rose"
#define ROSE_OUT "rose_outfile.txt"
#define ROSE_IN "rose_infile.txt"

#define GNUPLOT "~/10giga_volume/gnuplot/bin/gnuplot"
using namespace std;

void
createRandomTreeUsingBeep(SequenceTree &tree, int numLeafs, int ultrametricDeviation){
  //*********************
  //RUN beep in file
  remove(BEEP_OUT);
  
  //execute beep
  //creates an ultrametric birth death tree with hight 1.0.
  //  if( ! file_exists(BEEP) ){
  //  USER_WARNING("BEEP doesn't exist on the expected location: " BEEP);
  //  return;
  //}
  string beep_str = BEEP;
  beep_str += string(" -o " BEEP_OUT " -To params.true -Gn -Gt -Bp 0.4 0.4 ");
  beep_str = beep_str+numLeafs;
  
  cout<< " about to execute: " << beep_str<< endl;
  system(beep_str.c_str() );  
  //*********************
  //Read the tree string
  ifstream beep_outfile;
  open_read_stream(BEEP_OUT,beep_outfile);  
  tree = SequenceTree(beep_outfile);
  beep_outfile.close();

  SequenceTree::NodeVector nodes;
  tree.addNodesInPrefixOrder(nodes);
  EDGE(nodes[0]) = -1;//the root edge
  for(size_t i=1; i<nodes.size() ; i++){
    EDGE(nodes[i]) = EDGE(nodes[i])* randomFloat(((float)1)/ultrametricDeviation, 
						ultrametricDeviation);
  }
}

void
evolveSequencesUsingSeqGen(SequenceTree &tree, int seqlen, float diameterFactor, bool writeAncestral){

  // if( ! file_exists(SEQGEN) ){
//     USER_WARNING("SEQGEN doesn't exist on the expected location: " SEQGEN);
//     return;
//   }
  //Write to seq-gen file exucute SEQ-GEN
  ofstream seqgen_infile;
  remove(SEQGEN_IN);
  remove(SEQGEN_OUT);
  
  open_write_stream(SEQGEN_IN,seqgen_infile);
  seqgen_infile << tree << endl;
  seqgen_infile.close();

  //execute seq-gen
  //execute seq-gen
  char str[1000];
  sprintf(str,  SEQGEN " -l%d %s" //-wa =write ancestor sequences  
          //JUKES CANTOR
          //	  " -s %f -n 1 -m HKY -t 0.5 -f  0.25 0.25 0.25 0.25 < "
          //" -n1 -mHKY -t0.5 -f0.25,0.25,0.25,0.25 < " //without diamter factor
          //K2P with fix ratio 2
          " -s%f -n1 -mHKY -t2 -f0.25,0.25,0.25,0.25 < "
          SEQGEN_IN " > " SEQGEN_OUT,
          seqlen, (writeAncestral ? " -wa " : "") ,diameterFactor);
  cout << "about to execute: " << str << endl;
  system(str);

  ifstream seqgen_outfile;
  open_read_stream(SEQGEN_OUT,seqgen_outfile);
  

  if(!writeAncestral){
    tree.mapSequencesOntoTree(seqgen_outfile);
  }
  else{//seqgen writes the sequences in prefix order with left sibling first.
    vector<Sequence> seqs;    
    Sequence::readSequences(seqs,seqgen_outfile);
    SequenceTree::NodeVector nodes;
    tree.addNodesInPrefixOrder(nodes);
    ASSERT_EQ(seqs.size(),nodes.size());
    for(size_t i=0;i<nodes.size();i++){
      nodes[i]->data.s = seqs[i];
      if(nodes[i]->isLeaf()){
	ASSERT_EQ(NAME(nodes[i]),seqs[i].name);
      }
      else{
	NAME(nodes[i]) = string("i_") + NAME(nodes[i]);//better name of internal nodes
      }
    }
  }
  
  seqgen_outfile.close();
}

void
evolveSequencesUsingROSE(SequenceTree &tree,int seqlen, double exp_sub, double exp_indel){
  

//   if( ! file_exists(ROSE) ){
//     USER_WARNING("ROSE doesn't exist on the expected location: " ROSE);
//     return;
//  }
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
  tree.mapSequencesOntoTree(seqin);
  seqin.close();
}

void
createFileName(std::string &name, int numLeafs, int ultrametricDeviation, int seqlen, float diameterFactor){
  char tmp[1000];
  
  sprintf(tmp,"leafs_%d_sqlen_%d_diamF_%f.txt",
	  numLeafs,seqlen,diameterFactor);

  name.append(tmp);
}








void
createGnuplotFromDatFile(const char *datfilename, 
			 const string xlabel, const string ylabel,
			 const std::vector<int> columns, 
			 const std::vector<string> names, 
			 const std::string title){
  
  char plotfilename[1000];
  sprintf(plotfilename,"%s.eps", datfilename);

  
  ofstream out;
  system("rm -f ddd1234.txt");
  open_write_stream("ddd1234.txt",out);
  out.precision(2);
  out << "reset\n"
      <<"set terminal postscript eps solid color \"Time-Roman\" 18\n"
      <<"set output '" << plotfilename <<"'\n"
      <<"set xlabel '" << xlabel <<"'\n"
      <<"set ylabel '"<< ylabel<<"'\n"
      <<"set title '"<< title << "'\n"
      <<"set auto y\n"
      <<"set style data line\n"
      <<"set key left bottom\n"
      <<"set style fill solid border -1\n"
      <<"#set xtic rotate by 90\n"
      <<"#set bmargin 10\n";

  out <<"plot ";
  out << "'" << datfilename<<"' using 1:"<<columns[0] << " title '"<<names[0]<<"'";
  for(size_t i=1;i<columns.size();i++){
    out << ", '" << datfilename<<"' using 1:"<<columns[i] << " title '"<<names[i]<<"'";
  }
  out<<"\n"<< endl;
  out.close();
  string gnu(GNUPLOT " ddd1234.txt");

  cout<< " about to execute: " << gnu<< endl;
  system(gnu.c_str());
}



