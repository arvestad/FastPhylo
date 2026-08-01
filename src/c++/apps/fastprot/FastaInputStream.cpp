#include "FastaInputStream.hpp"

using namespace std;

FastaInputStream::FastaInputStream(char *filename)
    : reader(filename, "abcdefghiklmnopqrstuvwyzxABCDEFGHIKLMNOPQRSTUVWYZX -.?") {}

bool FastaInputStream::read(vector<Sequence> &seqs, string &runId, vector<string> &names, Extrainfos &extrainfos) {
  if (!readSequences(seqs, runId, extrainfos)) {
    return false;
}
  names.clear();
  names.reserve(seqs.size());
  for (size_t i = 0; i < seqs.size(); i++) {
    names.push_back(seqs[i].name);
  }
  return true;
}

bool FastaInputStream::readSequences(vector<Sequence> &seqs, string &runId, Extrainfos &extrainfos) {
  return reader.readSequences(seqs, runId, extrainfos);
}
