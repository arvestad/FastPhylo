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

#include "fastphylo/core/Exception.hpp"

#include <utility>
using namespace std;

Exception::Exception(std::string f, std::string func, int l, std::string mes) :
  file(std::move(f)),function(std::move(func)),line(l),message(std::move(mes)){
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





