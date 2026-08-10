#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

// SIMD-accelerated comparison primitives operating on ProtSeqCode-encoded
// sequences (see ProtSeqCode.hpp). Drop-in replacements for
// count_id_dist()/count_replacements() (ProtSeqUtils.cpp), operating on
// uint8_t residue codes instead of raw characters + toupper(). Not wired
// into the distance computation yet - that's Phase 3.
//
// Length-mismatch handling: count_id_dist() and count_replacements() both
// have undefined behavior if s2 is shorter than s1 (their iterators walk
// past s2.end()). These functions instead compare up to min(len1,len2)
// and, for count_mismatches()/count_id_fraction() only, treat any excess
// length in s1 as mismatches (keeping the denominator == len1, matching
// count_id_dist()'s existing convention) - see phase1_design.md for the
// rationale. This turns UB into defined behavior; it is a new decision,
// not a preserved one.

namespace ProtSeqCode
{

//! Number of mismatching positions among the first min(len1,len2) codes,
//! plus (len1 - len2) if s2 is shorter than s1.
std::size_t count_mismatches(const std::uint8_t *s1, std::size_t len1, const std::uint8_t *s2, std::size_t len2);

//! Identity fraction matching count_id_dist()'s contract: matches / len1.
double count_id_fraction(const std::uint8_t *s1, std::size_t len1, const std::uint8_t *s2, std::size_t len2);

//! Tally of canonical-amino-acid-to-canonical-amino-acid replacements,
//! flattened row-major NUM_CANONICAL_AA x NUM_CANONICAL_AA (index
//! [a * NUM_CANONICAL_AA + b], matching count_replacements()'s
//! Matrix(20,20) indexing temp(a,b)). Only positions i < min(len1,len2)
//! where both codes are canonical amino acids (is_canonical_aa()) are
//! tallied - matches getAAInd()'s existing sentinel-100 exclusion.
std::vector<std::size_t> count_replacement_tally(const std::uint8_t *s1, std::size_t len1, const std::uint8_t *s2,
                                                 std::size_t len2);

} // namespace ProtSeqCode
