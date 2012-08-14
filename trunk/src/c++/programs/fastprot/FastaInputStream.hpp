#ifndef FASTAINPUTSTREAM_HPP
#define FASTAINPUTSTREAM_HPP

#include "DataInputStream.hpp"

#include <iostream>

class FastaInputStream : public DataInputStream
{
public:
  FastaInputStream(char * filename );
  ~FastaInputStream();

  virtual bool read( std::vector<Sequence> &seqs, std::string & runId, std::vector<std::string> &names, Extrainfos &extrainfos );
  virtual bool readSequences( std::vector<Sequence> &seqs, std::string & runId, Extrainfos &extrainfos );
protected:

  bool readSeq(std::vector<Sequence> &seqs, std::string &line, int linesRead);

  std::istream * fp;
  std::ifstream fin;
  bool file_was_opened;
};

#endif // FASTAINPUTSTREAM_HPP
