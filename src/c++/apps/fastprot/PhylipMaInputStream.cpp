#include "PhylipMaInputStream.hpp"
#include <cstdio>

using namespace std;

bool PhylipMaInputStream::read(vector<Sequence> &seqs, string &runId, vector<string> &names, Extrainfos &extrainfos) {
  Sequence::readSequences(seqs, *reader.fp);
  names.clear();
  names.reserve(seqs.size());
  for (size_t i = 0; i < seqs.size(); i++) {
    names.push_back(seqs[i].name);
  }
  return true;
}

bool PhylipMaInputStream::readSequences(vector<Sequence> &seqs, string &runId, Extrainfos &extrainfos) {
  return reader.readSequences(seqs, runId, extrainfos);
}
