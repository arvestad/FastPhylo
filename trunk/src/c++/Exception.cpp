//--------------------------------------------------
//                                        
// File: Exception.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: Exception.cpp,v 1.3 2006/12/26 15:18:48 isaac Exp $                                 
//
//--------------------------------------------------

#include "Exception.hpp"
using namespace std;

Exception::Exception(const std::string &f, const std::string &func, int l, const std::string &mes) :
  file(f),function(func),line(l),message(mes){
}
Exception::Exception(const Exception &exc) :
file(exc.file), function(exc.function), line(exc.line), message(exc.message), stackTrace(exc.stackTrace){
}

std::ostream&
Exception::printOn(std::ostream& os) const {
  os << stackTrace << endl;
  
  os << "-----------------\n";
  os << "Exception\n";
  os << "File:      " << file << "\n";
  os << "Function:  " << function << "\n";
  os << "Line:      " << line << "\n";
  os << "Message:   " << message << "\n";
  os << "-----------------" << endl;
  return os;
}





