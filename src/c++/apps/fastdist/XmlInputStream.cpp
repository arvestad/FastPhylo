#include "XmlInputStream.hpp"
#include <cstdio>

using namespace std;

bool XmlInputStream::read(vector<DNA_b128_String> &b128seqs, string &runId, vector<string> &names, Extrainfos &extrainfos) {
  vector<Sequence> seqs;
  if (!readSequences(seqs, runId, extrainfos)) {
    return false;
}
  names.clear();
  names.reserve(seqs.size());
  for (auto & seq : seqs) {
    names.push_back(seq.name);
  }
  Sequences2DNA_b128(seqs, b128seqs);
  return true;
}

bool XmlInputStream::readSequences(vector<Sequence> &seqs, string &runId, Extrainfos &extrainfos) {
  return reader.readSequences(seqs, runId, extrainfos);
}
