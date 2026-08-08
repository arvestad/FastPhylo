#pragma once

#include <cstdio>
#include "DataOutputStream.hpp"

class XmlOutputStream : public DataOutputStream
{
  public:
    XmlOutputStream(char *filename);
    ~XmlOutputStream() override;

  protected:
    void print(tree2int_map &tree2count, bool printCounts, std::string &runId, std::vector<std::string> &names,
               Extrainfos &extrainfos) override;
};
