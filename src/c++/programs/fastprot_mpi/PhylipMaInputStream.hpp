#ifndef PHYLIPMAINPUTSTREAM_HPP
#define PHYLIPMAINPUTSTREAM_HPP

#include "DataInputStream.hpp"

#include <iostream>
#include <fstream>

using namespace std;

class PhylipMaInputStream : public DataInputStream
{
public:
  PhylipMaInputStream(char * filename = NULL);
  ~PhylipMaInputStream();

  virtual bool read( std::vector<Sequence> &seqs, std::string & runId, std::vector<std::string> &names, Extrainfos &extrainfos );
  virtual bool readSequences( std::vector<Sequence> &seqs, std::string & runId, Extrainfos &extrainfos );

protected:
  istream * fp;
  ifstream fin;
  bool file_was_opened;
};

#endif // PHYLIPMAINPUTSTREAM_HPP
