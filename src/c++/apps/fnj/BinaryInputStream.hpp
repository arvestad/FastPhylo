/*
 * BinaryInputStream.hpp
 *
 * 	Created on: Dec 14, 2011
 *      Auther: Mehmood Alam Khan
 *       Email: malagori@kth.se
 */
#pragma once

#include "DataInputStream.hpp"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

class BinaryInputStream : public DataInputStream {
public:
  BinaryInputStream(char *filename);
  ~BinaryInputStream() override;
  // Not an override - see DataInputStream.hpp's 2016-06-14 comment;
  // only the StrDblMatrix overload below is part of the base interface.
  readstatus readDM(StrFloMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos);
  readstatus readDM(StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos) override;

protected:
  istream *fp;
  ifstream fin;
  bool file_was_opened;
  int newSize;
  bool input_was_read;
};

