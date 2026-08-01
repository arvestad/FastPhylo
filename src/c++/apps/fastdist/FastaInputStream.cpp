#include "FastaInputStream.hpp"

using namespace std;

FastaInputStream::FastaInputStream(char *filename)
    : reader(filename, "acgtumrwsykvhdbnxACGTUMRWSYKVHDBNX -.?") {}

bool FastaInputStream::read(vector<DNA_b128_String> &b128seqs, string &runId, vector<string> &names, Extrainfos &extrainfos) {
  vector<Sequence> seqs;
  if (!readSequences(seqs, runId, extrainfos)) {
    return false;
  }
  names.clear();
  names.reserve(seqs.size());
  for (size_t i = 0; i < seqs.size(); i++) {
    names.push_back(seqs[i].name);
  }
  Sequences2DNA_b128(seqs, b128seqs);
  return true;
}

bool FastaInputStream::readSequences(vector<Sequence> &seqs, string &runId, Extrainfos &extrainfos) {
  return reader.readSequences(seqs, runId, extrainfos);
}
