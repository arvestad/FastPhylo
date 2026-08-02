#pragma once

#include "fastphylo/io/DataOutputStream.hpp"

// Layout Phase C: shared between fastdist and fastprot. Before this,
// fastprot already had the batched, one-fwrite()-per-row implementation
// (printPHYLIPfastSD, from speed2026a's output-speedup work); fastdist
// was still calling a ~200-line hand-rolled printPHYLIPfast() global
// function that issued one fprintf()/fwrite() call per matrix entry.
// Both wrote byte-identical output (verified via RunExamples.sh's ex1/
// ex3), so this is the same content, just reaching the file in
// numNodes I/O calls instead of numNodes*(numNodes+1)/2.
class PhylipDmOutputStream : public DataOutputStream {
public:
  PhylipDmOutputStream(char *filename) : DataOutputStream(filename), headerWritten(false) {}
  void print(StrDblMatrix &dm) override;
  void printSD(StrDblMatrix &dm) override;
  void printStartRun(std::vector<std::string> &names, std::string &runId, Extrainfos &extrainfos) override { headerWritten = false; }
  // fastdist's bootstrap-streaming path writes one row at a time as it's
  // computed (mem_eff_flag), rather than building the whole matrix
  // first - fastprot never needed this (it always has the whole matrix
  // before printing). Not exercised by RunExamples.sh, so left with its
  // original per-entry formatting rather than risking an unverifiable
  // rewrite; still shared here rather than duplicated, since the
  // ONEDIGIT/TENDIGIT-based formatting it uses is identical to
  // printPHYLIPfastSD's.
  void printRow(StrFloRow &dm, std::string name, int row, bool mem_eff_flag) override;
  // fastdist's main.cpp calls printHeader() unconditionally before its
  // printRow() streaming loop, needing it to actually write the count
  // line (there's no other point at which that line would be written in
  // that call pattern). fastprot's main.cpp *also* calls printHeader()
  // unconditionally, before print() - but relies on it being a no-op,
  // since printPHYLIPfastSD() below already writes its own count line.
  // Both callers are real and neither is wrong; `headerWritten` (reset
  // in printStartRun(), which precedes both call patterns) lets each
  // call site get what it needs without the two colliding into a
  // doubled count line. Caught by RunExamples.sh's ex17 fixture, not by
  // reasoning about it in advance.
  void printHeader(size_t numNodes) override;
  // printBootstrapSpliter: fastdist's original body was a single
  // commented-out line - already a no-op identical to the base class
  // default, so not overridden here.

  // skipLeadingHeader: when the plain-phylip "%5lu\n" node-count line has
  // already been written separately (by printHeader(), see above), skip
  // writing it again here. Only meaningful for the non-XML case; XmlOutputStream
  // never sets it (its printHeader() stays a no-op, so there's no
  // "already written" scenario to avoid there).
  enum class Format { Plain, Xml, XmlSD };
  static void printPHYLIPfastSD(const StrDblMatrix &dm, FILE *out, Format format, bool skipLeadingHeader = false);

private:
  bool headerWritten;
};
