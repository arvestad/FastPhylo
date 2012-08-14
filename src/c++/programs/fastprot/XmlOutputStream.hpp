#ifndef XMLOUTPUTSTREAM_HPP
#define XMLOUTPUTSTREAM_HPP

#include <cstdio>
#include "DataOutputStream.hpp"

class XmlOutputStream : public DataOutputStream
{
public:
  XmlOutputStream();
  XmlOutputStream(char * filename );
  virtual ~XmlOutputStream();

  virtual void print( StrDblMatrix & dm );
  virtual void printSD( StrDblMatrix & dm );
  virtual void printStartRun( std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos );
  virtual void printEndRun();
};

#endif // XMLOUTPUTSTREAM_HPP
