#pragma once

#include "DataInputStream.hpp"
#include "fastphylo/io/FastaInputStream.hpp"

// Layout Phase C: parsing logic lives in the shared io::FastaSequenceReader
// (see include/fastphylo/io/FastaInputStream.hpp); this class only adapts
// it to fastprot's DataInputStream interface.
class FastaInputStream : public DataInputStream {
public:
  FastaInputStream(char *filename);

  bool read(std::vector<Sequence> &seqs, std::string &runId, std::vector<std::string> &names, Extrainfos &extrainfos) override;
  bool readSequences(std::vector<Sequence> &seqs, std::string &runId, Extrainfos &extrainfos) override;

private:
  FastaSequenceReader reader;
};
