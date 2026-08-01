// Exhaustive alphabet/encoding test for ProtSeqCode (Phase 1 of the
// speed2026a protein-distance speedup; see phase1_design.md).
//
// Not yet wired into the CMake build (see phase0_audit.md - code_tests/
// isn't built by CMakeLists.txt today); compile and run manually:
//   c++ -std=c++11 -I../programs/fastprot \
//     ProtSeqCode_test.cpp ../programs/fastprot/ProtSeqCode.cpp -o /tmp/t && /tmp/t

#undef NDEBUG
#include <cassert>
#include <cctype>
#include <iostream>
#include <set>
#include "fastphylo/protein/ProtSeqCode.hpp"

using namespace ProtSeqCode;

static void test_canonical_aa_roundtrip() {
  const std::string canonical = "ARNDCQEGHILKMFPSTWYV";
  for (std::size_t i = 0; i < canonical.size(); i++) {
    char c = canonical[i];
    std::uint8_t code = encode_residue(c);
    assert(code == i);
    assert(is_canonical_aa(code));
    assert(decode_residue(code) == c);
    // lower-case must map to the same code
    assert(encode_residue(static_cast<char>(std::tolower(c))) == code);
  }
}

static void test_non_canonical_letters_excluded_from_replacement_tally() {
  const std::string extra = "BOUXZ";
  for (char c : extra) {
    std::uint8_t code = encode_residue(c);
    assert(!is_canonical_aa(code));
    assert(encode_residue(static_cast<char>(std::tolower(c))) == code);
  }
}

static void test_punctuation() {
  assert(!is_canonical_aa(encode_residue('-')));
  assert(!is_canonical_aa(encode_residue(' ')));
  assert(!is_canonical_aa(encode_residue('.')));
  assert(!is_canonical_aa(encode_residue('?')));
  assert(decode_residue(encode_residue('-')) == '-');
  assert(decode_residue(encode_residue(' ')) == ' ');
  assert(decode_residue(encode_residue('.')) == '.');
  assert(decode_residue(encode_residue('?')) == '?');
}

static void test_all_codes_distinct_for_fasta_legal_alphabet() {
  // Every character FastaInputStream.cpp's validation regex allows
  // must map to a distinct code from every other one (case folded).
  const std::string fasta_legal =
    "abcdefghiklmnopqrstuvwyzx -.?";
  std::set<std::uint8_t> seen;
  std::set<char> seen_upper;
  for (char c : fasta_legal) {
    char up = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    std::uint8_t code = encode_residue(c);
    if (seen_upper.count(up) != 0u) { continue; // already checked this identity class
}
    seen_upper.insert(up);
    assert(seen.find(code) == seen.end()); // no collision
    seen.insert(code);
    assert(code != OTHER_CODE); // every FASTA-legal symbol has its own code
  }
  assert(seen.size() == 29); // 25 letters (a-z minus j) + 4 punctuation
}

static void test_other_bucket_for_out_of_alphabet_bytes() {
  // 'j' is explicitly excluded by FastaInputStream's validation regex.
  assert(encode_residue('J') == OTHER_CODE);
  assert(encode_residue('j') == OTHER_CODE);
  assert(encode_residue('#') == OTHER_CODE);
  assert(encode_residue('\0') == OTHER_CODE);
  assert(encode_residue('9') == OTHER_CODE);
}

static void test_encode_sequence_matches_per_char_encode() {
  std::string seq = "ARNDcqeg-hilKMFPstWYVXbz??..  JJ";
  std::vector<std::uint8_t> out;
  encode_sequence(seq, out);
  assert(out.size() == seq.size());
  for (std::size_t i = 0; i < seq.size(); i++) {
    assert(out[i] == encode_residue(seq[i]));
}
}

static void test_empty_sequence() {
  std::vector<std::uint8_t> out;
  encode_sequence("", out);
  assert(out.empty());
}

int main() { // NOLINT(bugprone-exception-escape) - a test binary crashing on an unexpected exception is a fine, loud failure mode.
  test_canonical_aa_roundtrip();
  test_non_canonical_letters_excluded_from_replacement_tally();
  test_punctuation();
  test_all_codes_distinct_for_fasta_legal_alphabet();
  test_other_bucket_for_out_of_alphabet_bytes();
  test_encode_sequence_matches_per_char_encode();
  test_empty_sequence();
  std::cout << "ProtSeqCode_test: all tests passed" << std::endl;
  return 0;
}
