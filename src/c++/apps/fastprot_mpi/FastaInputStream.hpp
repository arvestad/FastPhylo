#ifndef FASTAINPUTSTREAM_HPP
#define FASTAINPUTSTREAM_HPP

#include "DataInputStream.hpp"
#include <iostream>

using namespace std;

class FastaInputStream : public DataInputStream {
public:
  FastaInputStream(char * filename );
  ~FastaInputStream();
  bool read( std::vector<Sequence> &seqs, std::string & runId, std::vector<std::string> &names, Extrainfos &extrainfos );
  bool readSequences( std::vector<Sequence> &seqs, std::string & runId, Extrainfos &extrainfos );

protected:
  bool readSeq(std::vector<Sequence> &seqs, std::string &line, int linesRead);
  istream * fp;
  ifstream fin;
  bool file_was_opened;
};

#endif // FASTAINPUTSTREAM_HPP
