#include "fastphylo/io/FastaInputStream.hpp"
#include "fastphylo/core/Exception.hpp"
#include <cstdio>

using namespace std;

FastaSequenceReader::~FastaSequenceReader() {
  if (file_was_opened)
    fin.close();
}

FastaSequenceReader::FastaSequenceReader(char *filename, string allowedChars) : allowedChars(allowedChars) {
  file_was_opened = false;
  if (filename == nullptr)
    fp = &cin;
  else {
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

bool FastaSequenceReader::readSeq(vector<Sequence> &seqs, string &line, int linesRead) {
  linesRead++;
  seqs.resize(linesRead);
  string seqStr;
  Sequence &s = seqs[linesRead - 1];
  string::size_type findPos = line.find_first_of("\t ", 1);
  string::size_type nameEndPos;
  if (findPos == string::npos) {
    nameEndPos = line.size() - 1;
  } else {
    nameEndPos = findPos - 1;
  }
  s.name = line.substr(1, nameEndPos);

  bool readGreaterThan = false;
  while (!readGreaterThan && getline(*fp, line)) {
    if (line.size() > 0) {
      if (line[0] == '>') {
        readGreaterThan = true;
        readSeq(seqs, line, linesRead);
      } else {
        seqStr += line;
      }
    }
  }
  if (seqStr.size() == 0 || seqStr.find_first_not_of(allowedChars) != string::npos) {
    THROW_EXCEPTION("Malformed Fasta format\n");
    exit(EXIT_FAILURE);
  } else {
    Sequence &s = seqs[linesRead - 1];
    s.seq = seqStr;
  }
  return true;
}

bool FastaSequenceReader::readSequences(vector<Sequence> &seqs, string &runId, Extrainfos &extrainfos) {
  string line;
  while (getline(*fp, line)) {
    line.erase(line.find_last_not_of(" \n\r\t") + 1);
    if (line.size() > 0 && line[0] == '>') {
      readSeq(seqs, line, 0);
      return true;
    }
  }
  THROW_EXCEPTION("Malformed Fasta format\n");
  exit(EXIT_FAILURE);
}
