//--------------------------------------------------
//                                        
// File: Exception.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: Exception.hpp,v 1.4 2006/12/26 15:18:48 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef EXCEPTION_HPP
#define EXCEPTION_HPP

#include "Object.hpp"
#include <string>
#include <iostream>
#include <exception>
#include <sstream>

//
// Macros for catching and throwing
//

#define THROW_EXCEPTION(MES)  {std::ostringstream out; out << MES; throw Exception(__FILE__,__FUNCTION__,__LINE__,out.str());}
#define TRY_EXCEPTION() try{
#define CATCH_EXCEPTION() }catch(Exception exc){ std::cerr << exc <<std::endl;}
#define CATCH_RETHROW() }catch(Exception exc1){ Exception exc2(__FILE__,__FUNCTION__,__LINE__,""); exc2.addToStackTrace(exc1); throw exc2;}


class Exception : public Object , public std::exception
{
public:
  std::string file;
  std::string function;
  int line;
  
  std::string message;


std::string stackTrace;
  
Exception(const std::string &f, const std::string &func, int l, const std::string &mes);
Exception(const Exception &exc);

  void addToStackTrace(const Exception &exc){
    std::ostringstream out;
    exc.printOn(out);
    stackTrace += out.str();
  }
  virtual ~Exception() throw() {
  }

  virtual std::ostream& printOn(std::ostream& os) const;

};


#endif // EXCEPTION_HPP










