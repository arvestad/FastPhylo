/*
 * BinaryInputStream.hpp
 *
 * 	Created on: Dec 14, 2011
 *      Auther: Mehmood Alam Khan
 *       Email: malagori@kth.se
 */
#ifndef BINARYINPUTSTREAM_HPP
#define BINARYINPUTSTREAM_HPP

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
  ~BinaryInputStream();
  readstatus readDM(StrFloMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos);
  readstatus readDM(StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos);

protected:
  istream *fp;
  ifstream fin;
  bool file_was_opened;
  int newSize;
  bool input_was_read;
};

#endif // BINARYINPUTSTREAM_HPP
