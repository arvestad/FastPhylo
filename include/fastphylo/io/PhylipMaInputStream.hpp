#pragma once

#include "fastphylo/core/Sequence.hpp"
#include "fastphylo/io/Extrainfos.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// Layout Phase C (input side): shared PHYLIP multi-alignment parsing.
// readSequences() and the constructor were byte-identical between
// fastdist's and fastprot's copies (confirmed by diff before merging,
// not assumed). fastdist's read() (DNA_b128_String output) doesn't go
// through readSequences() at all - it calls the separate global
// DNA_b128_StringsFromPHYLIP() directly on the underlying stream, so
// `fp` stays accessible to it rather than private.
class PhylipSequenceReader {
public:
  PhylipSequenceReader(char *filename);
  ~PhylipSequenceReader();

  bool readSequences(std::vector<Sequence> &seqs, std::string &runId, Extrainfos &extrainfos);

  std::istream *fp;

protected:
  std::ifstream fin;
  bool file_was_opened;
};
