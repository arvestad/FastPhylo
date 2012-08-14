#ifndef XMLOUTPUTSTREAMM_HPP
#define XMLOUTPUTSTREAMM_HPP

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

#endif // XMLOUTPUTSTREAMM_HPP
