#pragma once

#include "DataInputStream.hpp"
#include "fastphylo/io/XmlInputStream.hpp"
#include "fastphylo/io/fileFormatSchema.hpp"

// Layout Phase C: parsing lives in the shared io::XmlSequenceReader (see
// include/fastphylo/io/XmlInputStream.hpp); this class only adapts it to
// fastprot's DataInputStream interface and supplies the protein sequence
// RelaxNG schema.
class XmlInputStream : public DataInputStream
{
  public:
    XmlInputStream(char *filename = nullptr) : reader(filename, fastphylo_prot_sequence_xml_relaxngstr)
    {
    }

    bool read(std::vector<Sequence> &seqs, std::string &runId, std::vector<std::string> &names,
              Extrainfos &extrainfos) override;
    bool readSequences(std::vector<Sequence> &seqs, std::string &runId, Extrainfos &extrainfos) override;

  private:
    XmlSequenceReader reader;
};
