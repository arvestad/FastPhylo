//--------------------------------------------------
//                                        
// File: file_utils.cpp
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: file_utils.cpp,v 1.4 2006/09/10 10:18:26 isaac Exp $                                 
//
//--------------------------------------------------

#include "file_utils.hpp"
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "log_utils.hpp"
#include "Exception.hpp"

using namespace std;

/* CHECK IF FILE EXISTS */
int
file_exists(const char *fname){
  FILE *ftmp = fopen(fname,"r");
  if ( ftmp != NULL ){
    fclose(ftmp);
    return 1;
  }
  else
    return 0;
}



FILE *
open_write_file_interactive(const char *fname){
  const char *mode ="w";
  if ( file_exists(fname) ){
    while ( true ){
      cout << "File exists: \"" << fname << "\"" << endl;
      cout << "Do you want to (w)rite or (a)ppend? "<<endl;
      string choice;
      cin >> choice;
      if ( choice == "w" ){
        mode = "w";
        break;
      }
      else if ( choice == "a" ){
        mode = "a";
        break;
      }
      cout << "bad choice" << endl;
    }
  }

  string tmp;
  FILE *ftmp = fopen(fname,mode);
  if ( ftmp == NULL ){
    cout << "Couldn't open file \"" << fname << "\"" << endl;
    cout << "What file do you want to write to? " << endl; 
    cin >> tmp;
    fname = tmp.c_str();
    return open_write_file(fname);
  }
  else
    return ftmp;
}



FILE *
open_write_file(const char *fname) throw(Exception){
  const char *mode ="w";

  FILE *ftmp = fopen(fname,mode);
  if ( ftmp == NULL ){
    THROW_EXCEPTION("Couldn't open file \"" << fname << "\"");
  }
  else
    return ftmp;
}

void
open_write_stream_interactive(const char *fname, ofstream &of){
  ofstream::openmode mode = ofstream::out; 
  if ( file_exists(fname) ){
    while ( true ){
      cout << "File exists: \"" << fname << "\"" << endl;
      cout << "Do you want to (w)rite or (a)ppend? "<<endl;
      string choice;
      cin >> choice;
      if ( choice == "w" ){
        mode = ofstream::out;;
        break;
      }
      else if ( choice == "a" ){
        mode = ofstream::app;
        break;
      }
      cout << "bad choice" << endl;
    }
  }
  
  of.open(fname,mode);
  
  if ( !of.good() ){
    of.close();
    of.clear();
    USER_WARNING("can't open file " << string(fname));
    return;
  }
}

void
open_write_stream(const string fname, ofstream &of) throw(Exception){
  open_write_stream(fname.c_str(),of);
}

void
open_write_stream(const char *fname, ofstream &of) throw(Exception){
  ofstream::openmode mode = ofstream::out; 
  
  of.open(fname,mode);
  
  if ( !of.good() ){
    of.close();
    of.clear();
    THROW_EXCEPTION("Can't open file " << fname);
  }
}


FILE *
open_read_file_interactive(const char *fname){
  string tmp;
  while ( true ){
    FILE *ftmp = fopen(fname,"r");
    if ( ftmp != NULL )
      return ftmp;

    cout << "File doesn't exist: \"" << fname << "\"" << endl;
    cout << "What file do you want to read? " << endl;
  
    cin >> tmp;
    fname = tmp.c_str();
  } 
}

FILE *
open_read_file(const char *fname) throw(Exception){
  FILE *ftmp = fopen(fname,"r");
  if ( ftmp != NULL )
    return ftmp;
  THROW_EXCEPTION("File doesn't exist: \"" << fname <<"\"");
}

void
open_read_stream_interactive(const char *fname, ifstream &fin){

  string tmp;
  while ( true ){
    fin.open(fname,ifstream::in);
    
    if ( fin.good() )
      return;
    fin.close();
    fin.clear();
    cout << "File doesn't exist: \"" << fname << "\"" << endl;
    cout << "What file do you want to read? " << endl;
  
    cin >> tmp;
    fname = tmp.c_str();
  } 
}


void
open_read_stream(const string fname, ofstream &fin) throw(Exception){
  open_read_stream(fname.c_str(),fin);
}

void
open_read_stream(const char *fname, ifstream &fin) throw(Exception){
  fin.open(fname,ifstream::in);
    
  if ( fin.good() )
    return;

  fin.close();
  fin.clear();

  THROW_EXCEPTION("File doesn't exist: \"" << fname << "\"");
}


void
skipWhiteSpace(FILE *f){
  char c = fgetc(f);
  while ( isspace(c) )
    c = fgetc(f);

  ungetc(c,f);
}

void
skipWhiteSpace(std::istream &in){
  in >> ws;
}


void
skipUntil(std::istream &in, char *chars){

  char c;
  in.get(c);
  while ( in.good() && strchr(chars,(char)c) == NULL ){
    in.get(c);
  }
  in.unget();
}

void
appendToken(std::istream &in, std::string &str){
  char c;
  in.get(c);
  while ( in.good() && !isspace(c) ){
    str.push_back(c);
    in.get(c);
  }
  in.unget();
}


void
appendUntil(std::istream &in, std::string &str,  char *chars){
  char c;
  in.get(c);
  while (  in.good() && strchr(chars,c) == NULL ){
    str.push_back(c);
    in.get(c);
  }
  in.unget();

  
}
