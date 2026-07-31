#ifndef _PROTSEQUTILS_HPP_
#define _PROTSEQUTILS_HPP_

#include <vector>

  // Forward Declarations
  class Matrix;
  class Sequence;

  //! Translates an amino acid to an index
  std::size_t getAAInd(char c);

  //! Removes all indels in the given sequences
  void remove_gaps(std::vector<Sequence> &sv);

  //! Counts all replacements from an amino acid to another
  Matrix count_replacements(const Sequence &s1, const Sequence &s2);

  // count_id_dist() (percentage identity) used to live here; replaced by
  // ProtSeqCode::count_id_fraction() (ProtSeqCode.hpp/ProtSeqCompare.hpp)
  // as of speed2026a Phase 6 - see phase0_audit.md/phase1_design.md.

  //! Performs bootstrapping
  void bootstrap_sequences(const std::vector<Sequence> &seqs, std::vector<Sequence> &bseqs); 
#endif
