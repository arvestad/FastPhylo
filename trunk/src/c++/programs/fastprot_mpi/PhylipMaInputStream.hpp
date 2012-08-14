#ifndef PHYLIPMAINPUTSTREAMM_HPP
#define PHYLIPMAINPUTSTREAMM_HPP

#include "DataInputStream.hpp"

#include <iostream>
#include <fstream>

class PhylipMaInputStream : public DataInputStream
{
public:
  PhylipMaInputStream(char * filename );
  ~PhylipMaInputStream();

  virtual bool read( std::vector<Sequence> &seqs, std::string & runId, std::vector<std::string> &names, Extrainfos &extrainfos );
  virtual bool readSequences( std::vector<Sequence> &seqs, std::string & runId, Extrainfos &extrainfos );

protected:
  std::istream * fp;
  std::ifstream fin;
  bool file_was_opened;
};

#endif // PHYLIPMAINPUTSTREAMM_HPP
