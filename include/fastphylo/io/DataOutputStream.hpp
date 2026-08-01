#pragma once

#include <array>
#include <cstdio>
#include <vector>
#include <string>
#include "fastphylo/core/Exception.hpp"
#include "fastphylo/core/DistanceMatrix.hpp"
#include "fastphylo/core/DistanceRow.hpp"
#include "fastphylo/io/Extrainfos.hpp"

// Layout Phase C (project_layout_plan.md): shared base for fastdist's and
// fastprot's output streams - was two independently-diverged copies
// (programs/fastdist/DataOutputStream.hpp, programs/fastprot/
// DataOutputStream.hpp). The virtual interface below is the union of
// both: fastprot never needed printRows()/printBootstrapSpliter() (it
// writes the whole matrix at once via print()); fastdist never needed
// printSD() (no standard-deviation output). Default-empty bodies (not
// pure virtual) so each concrete stream only overrides what it uses -
// matches fastprot's original style, which needed less per-subclass
// boilerplate than fastdist's all-pure-virtual style for the same
// interface shape.
//
// fnj's DataOutputStream is NOT merged here: it writes trees built from
// a distance matrix, not a distance matrix itself, and reading it
// alongside this one confirmed the difference is real, not naming drift
// (see project_layout_plan.md's Phase C audit).
class DataOutputStream {
public:
  DataOutputStream(char *filename);
  // fastdist's original DataOutputStream had no destructor at all (an
  // empty inline `{}`), so `fp` - opened via fopen() in the constructor
  // below whenever a filename is given - was never fclose()'d. A real
  // file-descriptor leak on every invocation that writes to a named
  // file, parallel to (but independent of) the BinaryDmOutputStream
  // leak fixed earlier in this branch (Phase 6). fastprot's version
  // already got this right; that's the version kept here.
  virtual ~DataOutputStream();

  virtual void print(StrDblMatrix &dm) {}
  virtual void printSD(StrDblMatrix &dm) {}
  virtual void printStartRun(std::vector<std::string> &names, std::string &runId, Extrainfos &extrainfos) {}
  virtual void printEndRun() {}
  virtual void printRows(const StrDblMatrix &dm) {}
  virtual void printRow(StrFloRow &dm, std::string name, int row, bool mem_eff_flag = false) {}
  virtual void printHeader(size_t numNodes) {}
  virtual void printBootstrapSpliter(size_t numNodes) {}

  // Fast float-to-PHYLIP-field lookup tables, shared by PhylipDmOutputStream
  // and XmlOutputStream's formatting code.
  static const std::array<char, 128> ONEDIGIT;
  static const std::array<char, 128> TENDIGIT;

protected:
  bool file_was_opened;
  FILE *fp;
};
