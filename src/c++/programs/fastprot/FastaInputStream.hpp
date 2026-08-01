#pragma once

#include "DataInputStream.hpp"
#include <iostream>

using namespace std;

class FastaInputStream : public DataInputStream {
public:
  FastaInputStream(char * filename );
  ~FastaInputStream() override;
  bool read( std::vector<Sequence> &seqs, std::string & runId, std::vector<std::string> &names, Extrainfos &extrainfos ) override;
  bool readSequences( std::vector<Sequence> &seqs, std::string & runId, Extrainfos &extrainfos ) override;

protected:
  bool readSeq(std::vector<Sequence> &seqs, std::string &line, int linesRead);
  istream * fp;
  ifstream fin;
  bool file_was_opened;
};

