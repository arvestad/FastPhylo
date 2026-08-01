#pragma once

#include "DataInputStream.hpp"
#include <iostream>
#include <fstream>

class PhylipMaInputStream : public DataInputStream
{
public:
  PhylipMaInputStream(char * filename = nullptr);
  ~PhylipMaInputStream() override;

  bool read( std::vector<DNA_b128_String> &b128_strings, std::string & runId, std::vector<std::string> &names, Extrainfos &extrainfos ) override;
  bool readSequences( std::vector<Sequence> &seqs, std::string & runId, Extrainfos &extrainfos ) override;

protected:
  std::istream * fp;
  std::ifstream fin;
  bool file_was_opened;
};

