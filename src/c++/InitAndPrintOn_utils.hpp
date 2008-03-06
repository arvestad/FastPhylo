//--------------------------------------------------
//                                        
// File: InitAndPrintOn_utils.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: InitAndPrintOn_utils.hpp,v 1.10 2006/12/19 09:24:36 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef INITANDPRINTON_UTILS_HPP
#define INITANDPRINTON_UTILS_HPP

#include <string>
#include <iostream>
#include "file_utils.hpp"
#include "stl_utils.hpp"
#include "Sequence.hpp"

//
// This file contains function objects for printing and reading 
// some common data types from streams.
//

//---------------------------------
// PRIMITIVE DATA TYPES
// 
// Works for all data types for which the "<< i" and "i <<" operators
// are defined.
//
// Data_init<int>
// Data_init<float>
// Data_init<double>
// Data_init<std::string>
// Data_init<Sequence_float>

//INIT
template <class DataType>
struct Data_init
{
  void operator()(std::istream &in, DataType &d) const{
    skipWhiteSpace(in);
    in >> d;
  }
};
template <class DataType>
struct empty_Data_init
{
  void operator()(std::istream &in, const DataType &d) const{
    PROG_ERROR("Empty init function called");
  }
};

// PRINT ON
template <class DataType>
struct Data_printOn{

  std::ostream& operator()(std::ostream &os, const DataType &i) const{
    os << i;
    return os;
  }
};
template <class DataType>
struct empty_Data_printOn{

  std::ostream& operator()(std::ostream &os, const DataType &i) const{
    return os;
  }
};

//--------------------------------------
// SEQUENCE DOUBLE PAIR
//
struct Sequence_double{
  Sequence s;
  double dbl;
};


std::istream &
operator>>(std::istream &in,Sequence_double &seqflt);

std::ostream&
operator<<(std::ostream & os, const Sequence_double &seqflt);

//--------------------------------------
// STRING INT PAIR
//
struct string_int{
  std::string s;
  int i;
};

std::istream &
operator>>(std::istream &in,string_int &strint);

std::ostream&
operator<<(std::ostream & os, const string_int &strint);

//--------------------------------------
// STRING DOUBLE PAIR
//
struct string_double{
  std::string s;
  double dbl;
};


std::istream &
operator>>(std::istream &in,string_double &strint);

std::ostream&
operator<<(std::ostream & os, const string_double &strint);


//--------------------------------------
// INT DOUBLE PAIR
//
struct int_double{
  int i;
  double dbl;
};


std::istream &
operator>>(std::istream &in,int_double &intdbl);

std::ostream&
operator<<(std::ostream & os, const int_double &intdble);

#endif // INITANDPRINTON_UTILS_HPP















