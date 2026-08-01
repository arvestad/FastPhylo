#include "XmlInputStream.hpp"
#include <cstdio>

using namespace std;

bool XmlInputStream::read(vector<Sequence> &seqs, string &runId, vector<string> &names, Extrainfos &extrainfos) {
  if (!readSequences(seqs, runId, extrainfos)) {
    return false;
}
  names.clear();
  names.reserve(seqs.size());
  for (auto & seq : seqs) {
    names.push_back(seq.name);
  }
  return true;
}

bool XmlInputStream::readSequences(vector<Sequence> &seqs, string &runId, Extrainfos &extrainfos) {
  return reader.readSequences(seqs, runId, extrainfos);
}
