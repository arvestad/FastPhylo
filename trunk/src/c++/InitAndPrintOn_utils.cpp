//--------------------------------------------------
//                                        
// File: InitAndPrintOn_utils.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: InitAndPrintOn_utils.cpp,v 1.6 2006/12/08 11:09:12 isaac Exp $                                 
//
//--------------------------------------------------

#include "InitAndPrintOn_utils.hpp"
#include <string>
#include "file_utils.hpp"
#include "Newick.hpp"





std::istream &
operator>>(std::istream &in,Sequence_double &strflt){
    //There are three cases:  1. "str:flt" 2. "str" 3. ":flt"
    skipWhiteSpace(in);

    //get the first part
    strflt.s = Sequence();
    while ( strchr(":),;",in.peek()) == NULL && in.peek() != EOF )
      strflt.s.name += in.get();

     
    //get the second part
    if ( in.peek() == ':' ){
      std::string fltstr;
      in.get();//skip :
      while ( strchr("), ;",in.peek()) == NULL && in.peek() != EOF )
        fltstr += in.get();
      strflt.dbl = atof(fltstr.c_str());
    }
    else {
      strflt.dbl = -1;
    }
    
    //if the last semicol of the tree
    if (in.peek() == ';' ){
      in.get();
    }
    return in;
}

std::ostream&
operator<<(std::ostream & os, const Sequence_double &strflt){
    strflt.s.printShort(os);
    if ( strflt.dbl != -1 ) os  << newickDelimiters.sequence_double_left.c_str()
         <<  strflt.dbl << newickDelimiters.sequence_double_right.c_str();
    return os;
}


std::istream &
operator>>(std::istream &in, string_double &strflt){
  //There are three cases:  1. "str:flt" 2. "str" 3. ":flt"
  skipWhiteSpace(in);

  //get the first part
  strflt.s = std::string();
  while ( strchr(":),;",in.peek()) == NULL && in.peek() != EOF )
    strflt.s += in.get();

     
  //get the second part
  if ( in.peek() == ':' ){
    std::string fltstr;
    in.get();//skip :
    while ( strchr("), ;",in.peek()) == NULL && in.peek() != EOF )
      fltstr += in.get();
    strflt.dbl = atof(fltstr.c_str());
  }
  else {
    strflt.dbl = -1;
  }
    
  //if the last semicol of the tree
  if (in.peek() == ';' ){
    in.get();
  }
  return in;
}

std::ostream&
operator<<(std::ostream & os, const string_double &strflt){
  os << strflt.s;
  if ( strflt.dbl != -1 ) os  << ":"<< strflt.dbl;
  return os;
}




std::istream &
operator>>(std::istream &in,int_double &intdbl){
  //There are three cases:  1. "str:flt" 2. "str" 3. ":flt"
  skipWhiteSpace(in);

  //get the first part
  std::string s;
  while ( strchr(":),;",in.peek()) == NULL && in.peek() != EOF )
    s += in.get();
  intdbl.i = atoi(s.c_str());
     
  //get the second part
  if ( in.peek() == ':' ){
    std::string fltstr;
    in.get();//skip :
    while ( strchr("), ;",in.peek()) == NULL && in.peek() != EOF )
      fltstr += in.get();
    intdbl.dbl = atof(fltstr.c_str());
  }
  else {
    intdbl.dbl = -1;
  }
    
  //if the last semicol of the tree
  if (in.peek() == ';' ){
    in.get();
  }
  return in;
}

std::ostream&
operator<<(std::ostream & os, const int_double &intdbl){
  os << intdbl.i;
  if ( intdbl.dbl != -1 ) os  << ":"<< intdbl.dbl;
  return os;
}

std::istream &
operator>>(std::istream &in,string_int &strint){
    //There are three cases:  1. "str:flt" 2. "str" 3. ":flt"
    skipWhiteSpace(in);

    //get the first part
    strint.s = "";
    while ( strchr(":),;",in.peek()) == NULL && in.peek() != EOF )
      strint.s += in.get();

     
    //get the second part
    if ( in.peek() == ':' ){
      std::string fltstr;
      in.get();//skip :
      while ( strchr("), ;",in.peek()) == NULL && in.peek() != EOF )
        fltstr += in.get();
      strint.i = atoi(fltstr.c_str());
    }
    else {
      strint.i = -1;
    }
    
    //if the last semicol of the tree
    if (in.peek() == ';' ){
      in.get();
    }
    return in;
}

std::ostream&
operator<<(std::ostream & os, const string_int &strint){
  os << strint.s;
  if ( strint.i != -1 ) os  << ":"<< strint.i;
  return os;
}














