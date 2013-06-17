#ifndef DATAOUTPUTSTREAM_HPP
#define DATAOUTPUTSTREAM_HPP

#include "Exception.hpp"
#include "SequenceTree.hpp"
#include "Extrainfos.hpp"

#include <iostream>
#include <vector>
#include <string>

using namespace std;

class DataOutputStream {
public:
  DataOutputStream(char *filename);
  virtual ~DataOutputStream() {};
  virtual void print(tree2int_map & tree2count, bool printCounts, string & runId,vector<string> & names, Extrainfos & extrainfos) {};
protected:
  ostream * fp;
  ofstream fout;
  bool file_was_opened;
};

#endif // DATAOUTPUTSTREAM_HPP
