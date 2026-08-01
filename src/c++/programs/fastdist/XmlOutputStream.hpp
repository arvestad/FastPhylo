#pragma once

#include <cstdio>
#include "DataOutputStream.hpp"

class XmlOutputStream : public DataOutputStream
{
public:
  XmlOutputStream();
  XmlOutputStream(char * filename = nullptr);
  ~XmlOutputStream() override;

  void print( StrDblMatrix & dm ) override;
  void printStartRun( std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos ) override;
  void printEndRun() override;
  //mehmood changes here... email: malagori@kth.se
  void printRow( StrFloRow & dm , std::string name, int row, bool mem_eff_flag) override;
  void printHeader( size_t numNodes ) override;
  void printBootstrapSpliter(size_t numNodes) override;
};

