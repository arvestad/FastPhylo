#pragma once

#include <cstdio>
#include "DataOutputStream.hpp"
#include "PhylipDmOutputStream.hpp"

class XmlOutputStream : public DataOutputStream
{
public:
  XmlOutputStream();
  XmlOutputStream(char * filename = nullptr);
  ~XmlOutputStream() override;

  void print( StrDblMatrix & dm ) override;
  void printSD( StrDblMatrix & dm ) override;
  void printStartRun( std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos ) override;
  void printEndRun() override;
};

