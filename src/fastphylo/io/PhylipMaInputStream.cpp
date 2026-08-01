#include "fastphylo/io/PhylipMaInputStream.hpp"
#include "fastphylo/core/Exception.hpp"
#include <cstdio>

using namespace std;

PhylipSequenceReader::~PhylipSequenceReader() {
  if (file_was_opened) {
    fin.close();
}
}

PhylipSequenceReader::PhylipSequenceReader(char *filename) {
  file_was_opened = false;
  if (filename == nullptr) {
    fp = &std::cin;
  } else {
    fin.open(filename, ifstream::in);
    if (!fin.good()) {
      fin.close();
      fin.clear();
      THROW_EXCEPTION("File doesn't exist: \"" << filename << "\"");
    }
    file_was_opened = true;
    fp = &fin;
  }
}

bool PhylipSequenceReader::readSequences(std::vector<Sequence> &seqs, std::string &runId, Extrainfos &extrainfos) {
  Sequence::readSequences(seqs, *fp);
  return true;
}
