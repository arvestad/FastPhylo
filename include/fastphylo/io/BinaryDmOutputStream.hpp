#pragma once

#include "fastphylo/io/DataOutputStream.hpp"
#include <cstdio>
#include <fstream>
#include <memory>

// Layout Phase C: shared between fastdist and fastprot. Both had already
// independently fixed the same class of bug in this constructor
// (fastprot: Phase 2, `delete fp` on a still-open, never-nulled pointer -
// heap corruption; fastdist: Phase 6, `fp = nullptr; delete fp;` - a
// file-descriptor leak instead, since the null check made the delete
// harmless but the file was never closed). Both fixes converged on the
// same unique_ptr<ofstream>-based design kept here.
//
// print() (fastprot: writes the whole matrix at once) and printRow()
// (fastdist: streams one row at a time as it's computed, for the
// bootstrap path) are both kept - genuinely different use cases, not
// duplicate implementations of the same thing.
class BinaryDmOutputStream : public DataOutputStream {
public:
  // matricesPerRun: how many print()/printRow() matrix bodies follow
  // each printHeader() call before the next run's header (or EOF) -
  // binary_dm_format_plan.md's run-boundary fix. Both real callers
  // (fastprot/fastdist) already know this as a single value fixed for
  // the whole program invocation: (no_incl_orig ? 0 : 1) + bootstraps.
  BinaryDmOutputStream(char *filename, size_t matricesPerRun);
  void printHeader(size_t numNodes) override;
  void printStartRun(std::vector<std::string> &names, std::string &runId, Extrainfos &extrainfos) override;
  void print(StrDblMatrix &dm) override;
  void printRow(StrFloRow &dm, std::string name, size_t row, bool mem_eff_flag) override;

private:
  // Owns the file when writing to disk; unused (writeToCout) when
  // writing to stdout, in which case ofs aliases &cout instead.
  std::unique_ptr<std::ofstream> file_stream;
  std::ostream *ofs;
  std::vector<std::string> m_names;
  size_t m_matricesPerRun;
  bool writeToCout;
};
