//--------------------------------------------------
//                                        
// File: Clustal2gaplessPhylip.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: Tamir2Phylip.cpp,v 1.2 2006/08/08 13:50:51 isaac Exp $                                 
//
//--------------------------------------------------



#include <string>
#include <time.h>
#include "arg_utils_ext.hpp"
#include "file_utils.hpp"
#include <iostream>
#include <vector>

#include "Sequence.cpp"
#include <fstream>
#include "stl_utils.hpp"

using namespace std; 

void
print_options(char *note = NULL){
  if ( note != NULL ){
    cout << "ERROR: " << note <<endl<< endl;
  }
  cout <<
    "OPTIONS           MAN/DEF   DESCRIPTION\n"
    " -o file                    Output for phylip file. Default is input file with the suffix .pylip\n"
    " -h or --help               Print this help message\n"
    " -V                         Print version\n";
  cout << endl;

  cout << "EXAMPLE USAGE" << endl;
  cout << " $./Clustal2gaplessPhylip clustal.aln" <<endl;
  cout << " $./Clustal2gaplessPhylip clustal.aln -o phylip.txt" <<endl;
  
  
  exit(1);
}




int
main( int argc, char ** argv){

  ifstream in;
  if  ( argc < 2 || HAS_OPTION("-h") || HAS_OPTION("--help") )
    print_options("no input file");
  open_read_stream(argv[1],in);

  ofstream out;
  string outname(argv[1]);
  outname += ".phylip";
  char *outfile_name = GET_OPTION_VAL("-o");
  if ( outfile_name != NULL )
    outname = string(outfile_name);
  open_write_stream(outname.c_str(),out);

  //----------------------------------

  str2str_hashmap name2seq(100);
  vector<string> names;

  //skip the first line
  skipUntil(in,"\n");
  in.get();
  
  //read lines and skip commented lines and empty lines.
  while ( in.good() ){
    char c = in.peek();
    if ( isspace(c) || c == '#' || in.peek() == EOF){
      LINE();
      skipUntil(in,"\n");
      in.get();
      continue;
    }
    LINE();
    string name;
    appendToken(in,name);
    PRINT(name);
    skipWhiteSpace(in);
    //look up in hash
    str2str_hashmap::iterator iter = name2seq.find(name);
    if ( iter == name2seq.end() ){
      names.push_back(name);
      name2seq[name] = "";
      iter = name2seq.find(name);
    }
    string &seq = (*iter).second;
    in.get(c);
    while ( in.good() && c != '\n' ){
      if ( ! isspace(c) )
        seq.push_back(c);
      in.get(c);
    }
    in.unget();
    //    appendToken(in,seq);

    //skip until new line
    skipUntil(in,"\n");
  }

  //FINNISHED READING
  LINE();
  //check lengths
  vector<string> seqs;
  size_t len = 0;
  for ( size_t i = 0 ; i < names.size() ; i++ ){
        LINE();

    seqs.push_back((*name2seq.find(names[i])).second);
    if ( len != 0 ){
      if ( seqs.back().length() != len ){
        PRINT(names[i]);
        USER_ERROR("length of sequences not equal " << len << " != " << seqs.back().length());
      }
    }
    len = seqs.back().length();
  }
  
  //find columns with gaps.
  vector<bool> gaps(len,false);
  for ( size_t i = 0 ; i < seqs.size() ; i++ ){
    LINE();
    string &seq = seqs[i];
    for ( size_t j  = 0 ; j < len ; j++ ){
      if ( seq[j] == '-' )
        gaps[j] = true;
    }
  }
  LINE();
  //remove columns with gaps
  for ( size_t i = 0 ; i < seqs.size() ; i++ ){
    string &oldseq = seqs[i];
    string newseq;
    newseq.reserve(len);
    for ( size_t j  = 0 ; j < len ; j++ ){
      if ( ! gaps[j] )
        newseq.push_back(oldseq[j]);
    }
    seqs[i] = newseq;
  }
  //----------------------------

  //WRITE OUTPUT
  out << "   " << seqs.size() << "  " << seqs[0].length() << endl;
  for ( size_t i = 0 ; i < seqs.size() ; i++ ){
    if ( names[i].length() < 10 ){
      out << std::setw(10) << std::left;
      out << names[i];
    }
    else
      out << names[i] << " ";
    out << seqs[i] << endl;
  }
  

  //----------------------------
  in.close();
  out.close();
  return 1;
}
