#pragma once

#include "DataInputStream.hpp"
#include "fastphylo/io/FastaInputStream.hpp"

// Layout Phase C: parsing logic lives in the shared io::FastaSequenceReader
// (see include/fastphylo/io/FastaInputStream.hpp); this class only adapts
// it to fastdist's DataInputStream interface (read() returning
// DNA_b128_String, which fastprot's equivalent class doesn't need).
class FastaInputStream : public DataInputStream {
public:
  FastaInputStream(char *filename = nullptr);

  bool read(std::vector<DNA_b128_String> &b128seqs, std::string &runId, std::vector<std::string> &names, Extrainfos &extrainfos) override;
  bool readSequences(std::vector<Sequence> &seqs, std::string &runId, Extrainfos &extrainfos) override;

private:
  FastaSequenceReader reader;
};
