#ifndef DATAINPUTSTREAM_HPP
#define DATAINPUTSTREAM_HPP

#include <cstdio>
#include "Sequences2DistanceMatrix.hpp"
#include "Exception.hpp"

#include <iostream>
#include <fstream>

typedef enum { DM_READ = 1, END_OF_RUN = 2, END_OF_RUNS = 3, ERROR = 4 } readstatus;

class DataInputStream
{
public:
  DataInputStream();
  virtual ~DataInputStream() {};


  virtual readstatus readDM( StrDblMatrix & dm ) = 0;

};

class PhylipMaInputStream : public DataInputStream
{
public:
  PhylipMaInputStream(char * filename );
  ~PhylipMaInputStream();

  virtual readstatus readDM( StrDblMatrix & dm );

protected:
  std::istream * fp;
  std::ifstream fin;
  bool file_was_opened;

 
};


#endif // DATAINPUTSTREAM_HPP
