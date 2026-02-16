#ifndef PHYLIPDMINPUTSTREAM_HPP
#define PHYLIPDMINPUTSTREAM_HPP
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
  ~PhylipDmInputStream();
  readstatus readDM( StrDblMatrix & dm, vector<string> & names, string & runId, Extrainfos & extrainfos );
protected:
  istream * fp;
  ifstream fin;
  bool file_was_opened;
};

#endif // PHYLIPDMINPUTSTREAM_HPP
