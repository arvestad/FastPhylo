#include "XmlInputStream.hpp"
#include <cstdio>

using namespace std;

bool XmlInputStream::read(vector<Sequence> &seqs, string &runId, vector<string> &names, Extrainfos &extrainfos) {
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

bool XmlInputStream::readSequences(vector<Sequence> &seqs, string &runId, Extrainfos &extrainfos) {
  return reader.readSequences(seqs, runId, extrainfos);
}
