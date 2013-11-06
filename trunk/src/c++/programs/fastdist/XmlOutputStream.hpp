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
  virtual void printStartRun( std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos );
  virtual void printEndRun();
  //mehmood changes here... email: malagori@kth.se
  virtual void printRow( StrFloRow & dm , std::string name, int row, bool mem_eff_flag);
  virtual void printHeader( size_t numNodes );
  virtual void printBootstrapSpliter(size_t numNodes);
};

#endif // XMLOUTPUTSTREAM_HPP
