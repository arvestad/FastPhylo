#pragma once
/*
 * PhylipDmInputStream.hpp
 *      Auther: Mehmood Alam Khan Email: malagori@kth.se
 */
#include <cstdio>
#include "DataInputStream.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

class PhylipDmInputStream : public DataInputStream {
public:
  PhylipDmInputStream(char * filename);
  ~PhylipDmInputStream() override;
  readstatus readDM( StrDblMatrix & dm, vector<string> & names, string & runId, Extrainfos & extrainfos ) override;
protected:
  istream * fp;
  ifstream fin;
  bool file_was_opened;
};

