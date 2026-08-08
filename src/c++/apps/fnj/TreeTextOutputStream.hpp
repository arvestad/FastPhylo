#pragma once

#include "DataOutputStream.hpp"
#include "fastphylo/core/Exception.hpp"
#include "fastphylo/core/SequenceTree.hpp"
#include "fastphylo/io/Extrainfos.hpp"

#include <iostream>
#include <vector>
#include <string>

class TreeTextOutputStream : public DataOutputStream
{
  public:
    TreeTextOutputStream(char *filename) : DataOutputStream(filename) {};
    void print(tree2int_map &tree2count, bool printCounts, std::string &runId, std::vector<std::string> &names,
               Extrainfos &extrainfos) override;
};
