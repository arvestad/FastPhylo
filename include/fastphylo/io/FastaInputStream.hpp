#pragma once

#include "fastphylo/core/Sequence.hpp"
#include "fastphylo/io/Extrainfos.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// Layout Phase C (input side): shared FASTA-parsing implementation for
// fastdist and fastprot. Composed into (not inherited by) each program's
// own thin FastaInputStream, since the two domains' DataInputStream
// interfaces genuinely differ (fastdist's read() returns
// vector<DNA_b128_String>&, fastprot's returns vector<Sequence>& - see
// project_layout_plan.md's Phase C audit) and forcing that into one
// shared virtual base would need either templates or a second round of
// interface reconciliation for no real benefit; readSequences() (the
// part that's actually identical work) is what's shared here.
//
// fastdist's and fastprot's readSequences()/readSeq() had diverged into
// two structurally different algorithms (fastdist: one flat pass over
// every line, calling a separate readSeqName() on header lines and
// validating each body line incrementally as it's appended; fastprot:
// readSequences() finds the first header line and hands off to a
// recursive readSeq() that re-enters itself on each subsequent header,
// validating each sequence's whole accumulated body once). Kept
// fastprot's version (the one already exercised by more of this
// project's history) as the canonical implementation - verified
// byte-identical to fastdist's original output on both a single-line-
// per-sequence fixture (RunExamples.sh's ex2) and a hand-constructed
// multi-line-per-sequence file (not covered by RunExamples.sh, checked
// separately since it's exactly where the two algorithms' structural
// difference could have shown up).
class FastaSequenceReader {
public:
  // allowedChars: the residue alphabet accepted in a sequence body
  // (differs between DNA and protein - passed in rather than
  // hardcoded).
  FastaSequenceReader(char *filename, std::string allowedChars);
  ~FastaSequenceReader();

  bool readSequences(std::vector<Sequence> &seqs, std::string &runId, Extrainfos &extrainfos);

protected:
  bool readSeq(std::vector<Sequence> &seqs, std::string &line, int linesRead);

  std::istream *fp;
  std::ifstream fin;
  bool file_was_opened;
  std::string allowedChars;
};
