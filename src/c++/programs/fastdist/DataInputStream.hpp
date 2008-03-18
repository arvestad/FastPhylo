#ifndef DATAINPUTSTREAM_HPP
#define DATAINPUTSTREAM_HPP

#include <cstdio>
#include "Sequences2DistanceMatrix.hpp"
#include "Exception.hpp"
#include "Sequence.hpp"

#include <iostream>
#include <fstream>

class DataInputStream
{
public:
  DataInputStream();
  virtual ~DataInputStream() {};
  virtual bool read( std::vector<std::string> &names, std::vector<DNA_b128_String> &b128_strings ) = 0;
  virtual bool readSequences( std::vector<Sequence> &seqs ) = 0;
};

class PhylipMaInputStream : public DataInputStream
{
public:
  PhylipMaInputStream(char * filename );
  ~PhylipMaInputStream();

  bool read( std::vector<std::string> &names, std::vector<DNA_b128_String> &b128_strings);
  bool readSequences(std::vector<Sequence> &seqs);
protected:
  std::istream * fp;
  std::ifstream fin;
  bool file_was_opened;

};


#endif // DATAINPUTSTREAM_HPP
