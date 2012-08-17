/*
 * BinaryInputStream.hpp
 *
 * 	Created on: Dec 14, 2011
 *      Auther: Mehmood Alam Khan
 *       Email: malagori@kth.se
 */
#ifndef BINARYINPUTSTREAM_HPP
#define BINARYINPUTSTREAM_HPP

#include <cstdio>
#include "DataInputStream.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>


class BinaryInputStream : public DataInputStream
{
public:
  BinaryInputStream(char * filename );
  virtual ~BinaryInputStream();
  virtual readstatus readDM( StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos ) {} ;
  virtual readstatus readFloatDM( StrFloMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos );

protected:
  std::istream * fp;
  std::ifstream fin;
  bool file_was_opened;

 
};


#endif // BINARYINPUTSTREAM_HPP
