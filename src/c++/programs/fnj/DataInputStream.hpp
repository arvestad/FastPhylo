#ifndef DATAINPUTSTREAM_HPP
#define DATAINPUTSTREAM_HPP

#include <cstdio>
#include "Sequences2DistanceMatrix.hpp"
#include "Exception.hpp"

#include <iostream>
#include <fstream>

class DataInputStream
{
public:
  DataInputStream();
  virtual ~DataInputStream() {};
  virtual bool read(  StrDblMatrix & dm,  str2int_hashmap & name2id  ) = 0;
};

class PhylipMaInputStream : public DataInputStream
{
public:
  PhylipMaInputStream(char * filename );
  ~PhylipMaInputStream();

  bool read(  StrDblMatrix & dm,  str2int_hashmap & name2id );
protected:
  std::istream * fp;
  std::ifstream fin;
  bool file_was_opened;

};


#endif // DATAINPUTSTREAM_HPP
