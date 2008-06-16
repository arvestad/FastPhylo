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

  virtual bool readSpeciesNamesAndDM( std::vector<std::string> & speciesnames, StrDblMatrix & dm ) = 0;
  virtual bool readDM( StrDblMatrix & dm ) = 0;

};

class PhylipMaInputStream : public DataInputStream
{
public:
  PhylipMaInputStream(char * filename );
  ~PhylipMaInputStream();

  virtual  bool readSpeciesNamesAndDM( std::vector<std::string> & speciesnames, StrDblMatrix & dm );
  virtual bool readDM( StrDblMatrix & dm );

protected:
  std::istream * fp;
  std::ifstream fin;
  bool file_was_opened;

};


#endif // DATAINPUTSTREAM_HPP
