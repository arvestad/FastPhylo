#ifndef PHYLIPMAINPUTSTREAM_HPP
#define PHYLIPMAINPUTSTREAM_HPP

#include "DataInputStream.hpp"

#include <iostream>
#include <fstream>

class PhylipMaInputStream : public DataInputStream
{
public:
  PhylipMaInputStream(char * filename );
  ~PhylipMaInputStream();

  virtual bool read( std::vector<std::string> &names, std::vector<DNA_b128_String> &b128_strings);
  virtual bool readSequences(std::vector<Sequence> &seqs);
protected:
  std::istream * fp;
  std::ifstream fin;
  bool file_was_opened;
};

#endif // PHYLIPMAINPUTSTREAM_HPP
