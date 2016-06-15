//--------------------------------------------------
//                                        
// File: log_utils.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: log_utils.hpp,v 1.12 2006/12/24 08:45:43 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef LOG_UTILS_HPP
#define LOG_UTILS_HPP

#include <iostream>
#include <assert.h>
#include <time.h>
// If the PRINT statements should be executed
#ifndef USE_PRINT
#define USE_PRINT 1
#endif

// PRINT
#ifndef PRINT
#define PRINT(EXP) \
if(USE_PRINT){ std::cout << #EXP << " = " << (EXP) << std::endl ;}
#endif

// PRINT_V
#ifndef PRINT_V
#define PRINT_V(EXP) \
if(USE_PRINT){ std::cout << __FILE__<< ":" << __LINE__ << "  (" <<  #EXP << ") = " << (EXP) << std::endl ;}
#endif

// PRINT_V
#ifndef PRINT_TIME
#define PRINT_TIME(EXP) \
if(USE_PRINT){ clock_t log_utils_time = clock();\
EXP;\
std::cout << __FILE__<< ":" << __LINE__ << "  (" <<  #EXP << ")  took time  " <<\
 double(clock()-log_utils_time)/(CLOCKS_PER_SEC/1000) << " ms"<< std::endl ;}
#endif

// PRINT_EXP only prints the expression
#ifndef PRINT_EXP
#define PRINT_EXP(EXP) \
if(USE_PRINT){  std::cout << "executing: " <<  #EXP  << std::endl; EXP;}
#endif

// LINE
#ifndef LINE
#define LINE() \
if(USE_PRINT){ std::cout << __FILE__<< ":" << __LINE__ << std::endl;}
#endif

// SEPARATOR
#ifndef SEPARATOR
#define SEPARATOR() \
if(USE_PRINT){ std::cout << __FILE__<< ":" << __LINE__ << "------------------------------" << std::endl;}
#endif

//ASSERT
#ifndef ASSERT
#ifndef NDEBUG
#define ASSERT(EXP, P1, P2)\
if( !(EXP) ){ std::cerr << "************\nASSERT FAILED\nfile: "<< __FILE__ <<"\nfunc: " << __FUNCTION__<< "\nline: " << __LINE__<< "\n"\
 << "Expression not true: " << #EXP << "\n" << #P1 << " == " << (P1) << "\n" << #P2 << " == " << (P2) <<"\n**********" << std::endl; assert(EXP);}
#else
#define ASSERT(EXP, P1, P2) 
#endif
#endif



#ifndef ASSERT_EQ
#ifndef NDEBUG
#define ASSERT_EQ(EXP1, EXP2)\
if( !((EXP1)==(EXP2)) ){ std::cerr << "************\nASSERT_EQ FAILED\nfile: "<< __FILE__ <<"\nfunc: " << __FUNCTION__<< "\nline: " << __LINE__<< "\n"\
 << #EXP1 << " != " << #EXP2<< "\n" << #EXP1 << " = " << (EXP1) << "\n" << #EXP2 << " = " << (EXP2) <<"\n**********" << std::endl; assert((EXP1)==(EXP2));}
#else
#define ASSERT_EQ(EXP1, EXP2) 
#endif
#endif

//ASSERT
#ifndef ASSERT_OP
#ifndef NDEBUG
#define ASSERT_OP(EXP1, OP,EXP2)\
if( !((EXP1)OP(EXP2)) ){ std::cerr << "************\nASSERT_OP FAILED\nfile: "<< __FILE__ <<"\nfunc: " << __FUNCTION__<< "\nline: " << __LINE__<< "\n"\
 << #EXP1 << " !" << #OP << " " << #EXP2<< "\n" << #EXP1 << " = " << (EXP1) << "\n" << #EXP2 << " = " << (EXP2) <<"\n**********" << std::endl; assert((EXP1)OP(EXP2));}
#else
#define ASSERT_EQ(EXP1, EXP2) 
#endif
#endif

// MEM_CHECK
#ifndef MEM_CHECK
#define MEM_CHECK(PTR) \
if((PTR)==NULL){ std::cerr << "************\nOUT OF MEMORY\nfile: "<< __FILE__ <<"\nfunc: " << __FUNCTION__<< "\nline: " << __LINE__<< "\n" << #PTR << " == NULL "<< "\n************"<<std::endl ;exit(1);}
#endif

// PROG_ERR
#ifndef PROG_ERROR
#define PROG_ERROR(EXP) \
if(true){ std::cerr << "************\nPROGRAM ERROR\nfile: "<< __FILE__ <<"\nfunc: " << __FUNCTION__<< "\nline: " << __LINE__<< "\n" << EXP << "\n************"<< std::endl ;exit(1);}
#endif

// USER_ERR
#ifndef USER_ERROR
#define USER_ERROR(EXP) \
  if(true){ std::cerr << "************\nUSER INPUT ERROR\nfile: "<< __FILE__ <<"\nfunc: " << __FUNCTION__<< "\nline: " << __LINE__<< "\n" << EXP << "\n************"<< std::endl ;std::exit(1);}
#endif

// USER_WARNING
#ifndef USER_WARNING
#define USER_WARNING(EXP) \
if(true){ std::cerr << "************\nUSER WARNING\nfile: "<< __FILE__ <<"\nfunc: " << __FUNCTION__<< "\nline: " << __LINE__<< "\n" << EXP << "\n************"<< std::endl ;}
#endif

// DEBUG
#ifndef USE_DEBUG
#define USE_DEBUG 0
#endif
#ifdef DEBUG
#warning "DEBUG defined else where / not compatiable"
#else
#define DEBUG(EXP) \
if(USE_DEBUG) std::cout << __FILE__ << ":"__LINE__ <<"  (" << #EXP <<") = " << EXP << std::endl;
#endif


#endif // LOG_UTILS_HPP








