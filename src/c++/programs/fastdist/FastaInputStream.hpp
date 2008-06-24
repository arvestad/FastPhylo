#ifndef FASTAINPUTSTREAM_HPP
#define FASTAINPUTSTREAM_HPP

#include "DataInputStream.hpp"

#include <iostream>
#include <fstream>

class FastaInputStream : public DataInputStream
{
public:
  FastaInputStream(char * filename );
  ~FastaInputStream();

  virtual bool read( std::vector<std::string> &names, std::vector<DNA_b128_String> &b128_strings);
  virtual bool readSequences(std::vector<Sequence> &seqs);
protected:

  bool readSeq(std::vector<Sequence> &seqs, std::string &line, int linesRead);

  std::istream * fp;
  std::ifstream fin;
  bool file_was_opened;
};

#endif // FASTAINPUTSTREAM_HPP
