#ifndef DATAOUTPUTSTREAMM_HPP
#define DATAOUTPUTSTREAMM_HPP

#include <cstdio>
#include <vector>
#include <string>
#include "../../Exception.hpp"
#include "../../DistanceMatrix.hpp"
#include "Extrainfos.hpp"


class DataOutputStream
{
public:
  DataOutputStream( );
  DataOutputStream(char * filename );
  virtual ~DataOutputStream() {};
  virtual void print( StrDblMatrix & dm ) = 0;
  virtual void printSD( StrDblMatrix & dm ) = 0;
  virtual void printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos ) {};
  virtual void printEndRun() {};

  //Mehmood's addition here
    // FAST PRINTING OF FLOATS
    static const char ONEDIGIT[128];
    static const char TENDIGIT[128];

protected:
  FILE * fp;
  bool file_was_opened;
};

void
printPHYLIPfastSD(const StrDblMatrix &dm, FILE *out, bool writeXml, bool writeXmlSD );

class PhylipDmOutputStream : public DataOutputStream
{
public:
  PhylipDmOutputStream(char * filename ) : DataOutputStream(filename) {};
  virtual ~PhylipDmOutputStream() {};
  virtual void print( StrDblMatrix & dm );
  virtual void printSD( StrDblMatrix & dm );
};

#endif // DATAOUTPUTSTREAMM_HPP
