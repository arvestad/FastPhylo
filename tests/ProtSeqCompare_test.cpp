// Correctness test for ProtSeqCompare's primitives (Phase 2 of the
// speed2026a protein-distance speedup).
//
// The oracle here is a standalone mirror of ProtSeqUtils.cpp's
// getAAInd()/count_id_dist()/count_replacements(), not the production
// functions themselves - kept dependency-free so this test doesn't need
// to link Sequence/Matrix and their transitive dependencies. Phase 4
// adds the real end-to-end differential test against the production
// code, wired in behind the old/new switch; this test's job is only to
// validate ProtSeqCompare's primitives in isolation against a
// known-correct reference.
//
// This test is wired into the CMake build as the ProtSeqCompare_test
// target. To compile and run manually instead:
//   c++ -std=c++17 -I../include ProtSeqCompare_test.cpp \
//     ../src/fastphylo/protein/ProtSeqCode.cpp \
//     ../src/fastphylo/protein/ProtSeqCompare.cpp -o /tmp/t && /tmp/t

#undef NDEBUG
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <string_view>
#include <vector>
#include "fastphylo/protein/ProtSeqCode.hpp"
#include "fastphylo/protein/ProtSeqCompare.hpp"

using namespace ProtSeqCode;

namespace {

constexpr std::string_view CANONICAL = "ARNDCQEGHILKMFPSTWYV";
// Every character FastaInputStream.cpp's validation regex allows, one of
// each case-fold identity class (see phase1_design.md).
constexpr std::string_view FULL_ALPHABET = "ARNDCQEGHILKMFPSTWYVBOUXZ-. ?";

std::size_t ref_getAAInd(char c) {
  std::string_view::size_type pos = CANONICAL.find(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
  return pos == std::string_view::npos ? 100 : pos;
}

// Mirrors count_id_dist() exactly, including its assumption that s2 is at
// least as long as s1 (real UB otherwise - only call this with len(s2) >=
// len(s1) in this test).
double ref_count_id_dist(const std::string &s1, const std::string &s2) {
  double id = 0;
  for (std::size_t i = 0; i < s1.size(); i++) {
    if (std::toupper(static_cast<unsigned char>(s1[i])) == std::toupper(static_cast<unsigned char>(s2[i]))) {
      id++;
}
}
  return id / static_cast<double>(s1.size());
}

// Mirrors count_replacements() exactly (same s2-not-shorter-than-s1
// assumption).
std::vector<std::size_t> ref_count_replacements(const std::string &s1, const std::string &s2) {
  std::vector<std::size_t> counts(std::size_t{20} * 20, 0);
  for (std::size_t i = 0; i < s1.size(); i++) {
    std::size_t c1 = ref_getAAInd(s1[i]);
    std::size_t c2 = ref_getAAInd(s2[i]);
    if (c1 != 100 && c2 != 100) {
      counts[(c1 * 20) + c2]++;
}
  }
  return counts;
}

void check_against_reference(const std::string &s1, const std::string &s2) {
  assert(s2.size() >= s1.size()); // reference's assumption

  std::vector<std::uint8_t> c1;
  std::vector<std::uint8_t> c2;
  encode_sequence(s1, c1);
  encode_sequence(s2, c2);

  double got_id = count_id_fraction(c1.data(), c1.size(), c2.data(), c2.size());
  double want_id = ref_count_id_dist(s1, s2);
  if (s1.empty()) {
    assert(std::isnan(got_id) && std::isnan(want_id));
  } else {
    assert(std::fabs(got_id - want_id) < 1e-12);
}

  std::vector<std::size_t> got_tally = count_replacement_tally(c1.data(), c1.size(), c2.data(), c2.size());
  std::vector<std::size_t> want_tally = ref_count_replacements(s1, s2);
  assert(got_tally == want_tally);
}

std::string random_seq(std::size_t len, std::mt19937 &rng) {
  std::string s;
  s.reserve(len);
  for (std::size_t i = 0; i < len; i++) {
    s += FULL_ALPHABET[rng() % FULL_ALPHABET.size()];
}
  return s;
}

// rand_r() is POSIX-only (no MSVC equivalent); std::mt19937 is standard
// C++11, portable, and reproducible given the same seed - a strict
// improvement for test code, not just a portability patch.
void test_random_pairs_various_lengths_and_mismatch_rates() {
  std::mt19937 rng(12345);
  const std::array<std::size_t, 11> lengths = { 0, 1, 2, 17, 31, 32, 33, 100, 347, 1000, 3001 };
  for (unsigned long len : lengths) {
    std::string s1 = random_seq(len, rng);
    // low mismatch: copy s1 then perturb a few positions
    std::string s2_low = s1;
    for (std::size_t k = 0; k < len / 10; k++) {
      s2_low[rng() % len] = FULL_ALPHABET[rng() % FULL_ALPHABET.size()];
}
    // high mismatch: independent random sequence
    std::string s2_high = random_seq(len, rng);
    check_against_reference(s1, s1); // identical
    if (len > 0) {
      check_against_reference(s1, s2_low);
      check_against_reference(s1, s2_high);
    } else {
      check_against_reference(s1, s2_low); // both empty
    }
  }
}

void test_every_alphabet_symbol_at_boundaries_and_runs() {
  for (char c : FULL_ALPHABET) {
    // symbol at first and last position of an otherwise-'A' sequence
    std::string s1 = "AAAAAAAAAA";
    s1[0] = c;
    std::string s2 = s1;
    s2[s1.size() - 1] = c;
    check_against_reference(s1, s2);
    // a run of the symbol
    std::string run(20, c);
    check_against_reference(run, run);
  }
}

void test_edge_cases() {
  check_against_reference("", "");         // empty
  check_against_reference("A", "A");       // single residue, identical
  check_against_reference("A", "R");       // single residue, different
  check_against_reference("AAAA", "AAAA"); // identical
  check_against_reference("AAAA", "RRRR"); // completely different
  check_against_reference("----", "----"); // all-gap identical
  check_against_reference("----", "AAAA"); // all-gap vs all-residue
}

void test_length_mismatch_defined_behavior() {
  // s2 shorter than s1: old code has UB here, so no reference call.
  // New behavior (phase1_design.md): compare up to min(len1,len2);
  // treat s1's excess tail as mismatches; denominator stays len1.
  std::vector<std::uint8_t> s1;
  std::vector<std::uint8_t> s2;
  encode_sequence("AAAAA", s1);  // len 5
  encode_sequence("AAA", s2);    // len 3, shares "AAA" prefix with s1

  std::size_t mismatches = count_mismatches(s1.data(), s1.size(), s2.data(), s2.size());
  assert(mismatches == 2); // positions 3,4 have no counterpart in s2

  double id = count_id_fraction(s1.data(), s1.size(), s2.data(), s2.size());
  assert(std::fabs(id - 3.0 / 5.0) < 1e-12); // 3 matches out of denominator 5

  // s2 longer than s1: matches count_id_dist()'s existing (defined, no UB)
  // behavior exactly - only the first len1 characters of s2 are compared.
  std::vector<std::uint8_t> s3;
  std::vector<std::uint8_t> s4;
  encode_sequence("AAA", s3);     // len 3
  encode_sequence("AAAAA", s4);   // len 5, first 3 chars match s3
  assert(count_mismatches(s3.data(), s3.size(), s4.data(), s4.size()) == 0);
}

} // anonymous namespace

int main() { // NOLINT(bugprone-exception-escape) - a test binary crashing on an unexpected exception is a fine, loud failure mode.
  test_random_pairs_various_lengths_and_mismatch_rates();
  test_every_alphabet_symbol_at_boundaries_and_runs();
  test_edge_cases();
  test_length_mismatch_defined_behavior();
  std::cout << "ProtSeqCompare_test: all tests passed" << std::endl;
  return 0;
}
