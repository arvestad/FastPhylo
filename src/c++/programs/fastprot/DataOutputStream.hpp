#ifndef DATAOUTPUTSTREAM_HPP
#define DATAOUTPUTSTREAM_HPP

#include <cstdio>
#include <vector>
#include <string>
#include <math.h>
#include "Exception.hpp"
#include "DistanceMatrix.hpp"
#include "DistanceRow.hpp"
#include "Extrainfos.hpp"

class DataOutputStream {
public:
  DataOutputStream( );
  DataOutputStream(char * filename );
  virtual ~DataOutputStream();
  virtual void print( StrDblMatrix & dm ) {};
  virtual void printSD( StrDblMatrix & dm ) {};
  virtual void printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos ) {};
  virtual void printEndRun() {};
  virtual void printRows(const StrDblMatrix &dm) {};
  virtual void printRow( StrFloRow & dm , std::string name, int row) {};
  virtual void printHeader( size_t numNodes ) {};

  //Mehmood's addition here
    // FAST PRINTING OF FLOATS
  static const char ONEDIGIT[128];
  static const char TENDIGIT[128];

protected:
  bool file_was_opened;
  FILE * fp;
};

#endif // DATAOUTPUTSTREAM_HPP
