#pragma once

#include <cstdio>
#include "DataOutputStream.hpp"
#include "PhylipDmOutputStream.hpp"

class XmlOutputStream : public DataOutputStream
{
public:
  XmlOutputStream();
  XmlOutputStream(char * filename = nullptr);
  ~XmlOutputStream();

  void print( StrDblMatrix & dm );
  void printSD( StrDblMatrix & dm );
  void printStartRun( std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos );
  void printEndRun();
};

