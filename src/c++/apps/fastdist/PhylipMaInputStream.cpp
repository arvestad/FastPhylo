#include "PhylipMaInputStream.hpp"
#include <cstdio>

using namespace std;

bool PhylipMaInputStream::read(vector<DNA_b128_String> &b128_strings, string &runId, vector<string> &names, Extrainfos &extrainfos) {
  DNA_b128_StringsFromPHYLIP(*reader.fp, names, b128_strings);
  return true;
}

bool PhylipMaInputStream::readSequences(vector<Sequence> &seqs, string &runId, Extrainfos &extrainfos) {
  return reader.readSequences(seqs, runId, extrainfos);
}
