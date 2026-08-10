#pragma once

#include "DataInputStream.hpp"
#include "fastphylo/io/PhylipMaInputStream.hpp"

// Layout Phase C: parsing lives in the shared io::PhylipSequenceReader
// (see include/fastphylo/io/PhylipMaInputStream.hpp); this class only
// adapts it to fastprot's DataInputStream interface.
class PhylipMaInputStream : public DataInputStream
{
  public:
    PhylipMaInputStream(char *filename = nullptr) : reader(filename)
    {
    }

    bool read(std::vector<Sequence> &seqs, std::string &runId, std::vector<std::string> &names,
              Extrainfos &extrainfos) override;
    bool readSequences(std::vector<Sequence> &seqs, std::string &runId, Extrainfos &extrainfos) override;

  private:
    PhylipSequenceReader reader;
};
