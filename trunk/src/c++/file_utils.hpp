//--------------------------------------------------
//                                        
// File: file_utils.hpp
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: file_utils.hpp,v 1.4 2006/09/10 10:18:26 isaac Exp $                                 
//
//--------------------------------------------------

#ifndef FILE_UTILS_HPP
#define FILE_UTILS_HPP


#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "Exception.hpp"
#include <string>

/* CHECK IF FILE EXISTS */
int
file_exists(const char *fname);


FILE *
open_write_file(const char *fname) throw(Exception);

void
open_write_stream(const char *fname, std::ofstream &fout) throw(Exception);

void
open_write_stream(const std::string fname, std::ofstream &fout) throw(Exception);

FILE *
open_write_file_interactive(const char *fname);

void
open_write_stream_interactive(const char *fname, std::ofstream &fout);

FILE *
open_read_file(const char *fname) throw(Exception);

void
open_read_stream(const char *fname,std::ifstream &fin) throw(Exception);

void
open_read_stream(const std::string fname,std::ifstream &fin) throw(Exception);


FILE *
open_read_file_interactive(const char *fname);

void
open_read_stream_interactive(const char *fname,std::ifstream &fin);

void
skipWhiteSpace(FILE *f);

void
skipWhiteSpace(std::istream &in);

void
skipUntil(std::istream &in, char *chars);

void
appendToken(std::istream &in, std::string &str);

void
appendUntil(std::istream &in, std::string &str,  char *chars);


#endif // FILE_UTILS_HPP
