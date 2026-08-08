#pragma once

#include <vector>

// Forward Declarations
class Sequence;

//! Removes all indels in the given sequences
void remove_gaps(std::vector<Sequence> &sv);

// getAAInd()/count_id_dist()/count_replacements() used to live here;
// replaced by ProtSeqCode::count_id_fraction()/count_replacement_tally()
// (ProtSeqCode.hpp/ProtSeqCompare.hpp) as of speed2026a Phase 6 (id) and
// the count_replacements wiring round (replacements) - see
// phase0_audit.md/phase1_design.md.

//! Performs bootstrapping
void bootstrap_sequences(const std::vector<Sequence> &seqs, std::vector<Sequence> &bseqs);
