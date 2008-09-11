#ifndef DATAOUTPUTSTREAM_HPP
#define DATAOUTPUTSTREAM_HPP

#include <cstdio>
#include <vector>
#include <string>
#include "Exception.hpp"
#include "DistanceMatrix.hpp"
#include "Extrainfos.hpp"


class DataOutputStream
{
public:
  DataOutputStream( );
  DataOutputStream(char * filename );
  virtual ~DataOutputStream() {};
  virtual void print( StrDblMatrix & dm ) = 0;
  virtual void printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos ) {};
  virtual void printEndRun() {};

protected:
  FILE * fp;
  bool file_was_opened;
};

class PhylipDmOutputStream : public DataOutputStream
{
public:
  PhylipDmOutputStream(char * filename ) : DataOutputStream(filename) {};
  virtual ~PhylipDmOutputStream() {};
  virtual void print( StrDblMatrix & dm );
};

#endif // DATAOUTPUTSTREAM_HPP
