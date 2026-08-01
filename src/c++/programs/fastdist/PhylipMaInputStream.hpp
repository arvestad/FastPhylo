#pragma once

#include "DataInputStream.hpp"
#include "fastphylo/io/PhylipMaInputStream.hpp"

// Layout Phase C: parsing lives in the shared io::PhylipSequenceReader
// (see include/fastphylo/io/PhylipMaInputStream.hpp); this class only
// adapts it to fastdist's DataInputStream interface. read() bypasses
// readSequences() entirely - it calls the DNA-specific
// DNA_b128_StringsFromPHYLIP() directly on the reader's stream.
class PhylipMaInputStream : public DataInputStream {
public:
  PhylipMaInputStream(char *filename = nullptr) : reader(filename) {}

  bool read(std::vector<DNA_b128_String> &b128_strings, std::string &runId, std::vector<std::string> &names, Extrainfos &extrainfos) override;
  bool readSequences(std::vector<Sequence> &seqs, std::string &runId, Extrainfos &extrainfos) override;

private:
  PhylipSequenceReader reader;
};
