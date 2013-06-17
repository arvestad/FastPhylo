#ifndef TREETEXTOUTPUTSTREAM_HPP
#define TREETEXTOUTPUTSTREAM_HPP

#include "DataOutputStream.hpp"
#include "Exception.hpp"
#include "SequenceTree.hpp"
#include "Extrainfos.hpp"

#include <iostream>
#include <vector>
#include <string>

class TreeTextOutputStream : public DataOutputStream {
public:
  TreeTextOutputStream(char *filename): DataOutputStream(filename) {};
  void print( tree2int_map & tree2count, bool printCounts, std::string & runId, std::vector<std::string> & names, Extrainfos & extrainfos);
};

#endif // TREETEXTOUTPUTSTREAM_HPP
