#ifndef DATAINPUTSTREAM_HPP
#define DATAINPUTSTREAM_HPP

#include <cstdio>
#include "Sequences2DistanceMatrix.hpp"
#include "Exception.hpp"
#include "Sequence.hpp"
#include "Extrainfos.hpp"


#include <iostream>
#include <fstream>
#include <libxml/tree.h>


class DataInputStream
{
public:
  DataInputStream() {};
  virtual ~DataInputStream() {};
  virtual bool read( std::vector<std::string> &names, Extrainfos &extrainfos, std::vector<DNA_b128_String> &b128_strings ) = 0;
  virtual bool readSequences( std::vector<Sequence> &seqs, Extrainfos &extrainfos ) = 0;
};

#endif // DATAINPUTSTREAM_HPP
