#include "ProtSeqCode.hpp"
#include <array>
#include <cctype>

namespace ProtSeqCode {

  namespace {

    std::array<std::uint8_t, 256> build_table() {
      std::array<std::uint8_t, 256> table;
      table.fill(OTHER_CODE);

      // Canonical amino acids, same order as getAAInd() (ProtSeqUtils.cpp).
      const char canonical[NUM_CANONICAL_AA] = {
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
        'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'
      };
      for (std::size_t i = 0; i < NUM_CANONICAL_AA; i++) {
        auto up = static_cast<unsigned char>(canonical[i]);
        auto lo = static_cast<unsigned char>(std::tolower(up));
        table[up] = static_cast<std::uint8_t>(i);
        table[lo] = static_cast<std::uint8_t>(i);
      }

      // Other FASTA-legal letters that aren't canonical amino acids.
      const char other_letters[] = { 'B', 'O', 'U', 'X', 'Z' };
      for (std::size_t i = 0; i < sizeof(other_letters); i++) {
        auto up = static_cast<unsigned char>(other_letters[i]);
        auto lo = static_cast<unsigned char>(std::tolower(up));
        auto code = static_cast<std::uint8_t>(NUM_CANONICAL_AA + i);
        table[up] = code;
        table[lo] = code;
      }

      table[static_cast<unsigned char>('-')] = 25;
      table[static_cast<unsigned char>(' ')] = 26;
      table[static_cast<unsigned char>('.')] = 27;
      table[static_cast<unsigned char>('?')] = 28;

      return table;
    }

    const std::array<std::uint8_t, 256> &code_table() {
      static const std::array<std::uint8_t, 256> table = build_table();
      return table;
    }

    // Canonical uppercase representative for each code, for decode_residue().
    std::array<char, ALPHABET_SIZE> build_decode_table() {
      std::array<char, ALPHABET_SIZE> table;
      const char symbols[ALPHABET_SIZE] = {
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
        'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
        'B', 'O', 'U', 'X', 'Z',
        '-', ' ', '.', '?',
        '?' // OTHER_CODE: no single canonical character, use '?'
      };
      for (std::size_t i = 0; i < ALPHABET_SIZE; i++)
        table[i] = symbols[i];
      return table;
    }

  } // anonymous namespace

  std::uint8_t encode_residue(char c) {
    return code_table()[static_cast<unsigned char>(c)];
  }

  char decode_residue(std::uint8_t code) {
    static const std::array<char, ALPHABET_SIZE> table = build_decode_table();
    if (code >= ALPHABET_SIZE)
      return '?';
    return table[code];
  }

  void encode_sequence(const std::string &seq, std::vector<std::uint8_t> &out) {
    out.resize(seq.size());
    const std::array<std::uint8_t, 256> &table = code_table();
    for (std::size_t i = 0; i < seq.size(); i++)
      out[i] = table[static_cast<unsigned char>(seq[i])];
  }

}
