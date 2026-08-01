#pragma once

#include "DataInputStream.hpp"
#include <iostream>
#include <fstream>

class PhylipMaInputStream : public DataInputStream
{
public:
  PhylipMaInputStream(char * filename = NULL);
  ~PhylipMaInputStream();

  virtual bool read( std::vector<DNA_b128_String> &b128_strings, std::string & runId, std::vector<std::string> &names, Extrainfos &extrainfos );
  virtual bool readSequences( std::vector<Sequence> &seqs, std::string & runId, Extrainfos &extrainfos );

protected:
  std::istream * fp;
  std::ifstream fin;
  bool file_was_opened;
};

