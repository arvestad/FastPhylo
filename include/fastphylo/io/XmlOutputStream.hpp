#pragma once

#include <cstdio>
#include "fastphylo/io/DataOutputStream.hpp"
#include "fastphylo/io/PhylipDmOutputStream.hpp"

// Layout Phase C: shared between fastdist and fastprot. print()/printSD()
// route through PhylipDmOutputStream::printPHYLIPfastSD() (writeXml=true)
// - fastprot already did this before this merge; fastdist's own print()
// called the old, unbatched printPHYLIPfast() global function instead.
class XmlOutputStream : public DataOutputStream {
public:
  XmlOutputStream(char *filename = nullptr);
  ~XmlOutputStream() override;

  void print(StrDblMatrix &dm) override;
  void printSD(StrDblMatrix &dm) override;
  void printStartRun(std::vector<std::string> &names, std::string &runId, Extrainfos &extrainfos) override;
  void printEndRun() override;

  // fastdist's memory-efficient bootstrap-streaming path: writes one row
  // at a time as it's computed, rather than building the whole matrix
  // first. fastprot never needed this. Not exercised by RunExamples.sh,
  // kept behaviorally identical to the original per-entry implementation
  // rather than risking an unverifiable rewrite (same reasoning as
  // PhylipDmOutputStream::printRow).
  void printRow(StrFloRow &dm, std::string name, int row, bool mem_eff_flag) override;
  void printBootstrapSpliter(size_t numNodes) override;
};
