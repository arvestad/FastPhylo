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
#pragma once

#include "fastphylo/core/Object.hpp"
#include <string>
#include <iostream>
#include <exception>
#include <sstream>
#include <cstdlib>

//
// Macros for catching and throwing
//

#define THROW_EXCEPTION(MES)  {std::ostringstream out; out << MES; throw Exception(__FILE__,__FUNCTION__,__LINE__,out.str());}
#define TRY_EXCEPTION() try{
// legacy_error_handling_plan.md Finding 0: an app that hit a real
// Exception used to still exit 0 (success) here - nothing after this
// catch block called exit(EXIT_FAILURE), so a caller checking the
// exit code was silently lied to.
#define CATCH_EXCEPTION() }catch(Exception exc){ std::cerr << exc <<std::endl; exit(EXIT_FAILURE);}
#define CATCH_RETHROW() }catch(Exception exc1){ Exception exc2(__FILE__,__FUNCTION__,__LINE__,""); exc2.addToStackTrace(exc1); throw exc2;}


class Exception : public Object , public std::exception
{
public:
  std::string file;
  std::string function;
  int line;
  
  std::string message;


std::string stackTrace;
  
Exception(std::string f, std::string func, int l, std::string mes);
Exception(const Exception &exc);

  void addToStackTrace(const Exception &exc){
    std::ostringstream out;
    exc.printOn(out);
    stackTrace += out.str();
  }
  // Modernization Phase 0 (modernization_plan.md): throw() (means
  // "never throws") is the pre-C++11 spelling of noexcept.
  ~Exception() noexcept override {
  }

  std::ostream& printOn(std::ostream& os) const override;

};












