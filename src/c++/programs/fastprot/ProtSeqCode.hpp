#ifndef _PROTSEQCODE_HPP_
#define _PROTSEQCODE_HPP_

#include <cstdint>
#include <string>
#include <vector>

// Byte-per-residue encoding of protein sequences, designed to let
// mismatch counting and replacement tallying operate on flat uint8_t
// arrays instead of re-deriving toupper()'d character identity on every
// comparison. See phase1_design.md for the full alphabet rationale.
//
// Codes 0-19: canonical amino acids, same order as getAAInd()
//             (ProtSeqUtils.cpp): A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V
// Codes 20-24: other FASTA-legal letters that aren't canonical AAs:
//              B,O,U,X,Z
// Code 25: '-' (gap/indel)
// Code 26: ' ' (space)
// Code 27: '.' (dot)
// Code 28: '?' (question mark)
// Code 29: OTHER - any byte not covered above (see phase1_design.md for
//          the disclosed behavior-change corner case this covers)

namespace ProtSeqCode {

  const std::size_t NUM_CANONICAL_AA = 20;
  const std::size_t ALPHABET_SIZE = 30;
  const std::uint8_t OTHER_CODE = 29;

  //! Encodes a single residue character (case-insensitive) to its code.
  std::uint8_t encode_residue(char c);

  //! Decodes a code back to its canonical uppercase character
  //! representative. Not performance-sensitive; for messages/debugging.
  char decode_residue(std::uint8_t code);

  //! Encodes a whole sequence, one code per residue, in order.
  void encode_sequence(const std::string &seq, std::vector<std::uint8_t> &out);

  //! True iff code identifies one of the 20 canonical amino acids
  //! (i.e. what getAAInd() would have accepted instead of sentinel 100).
  inline bool is_canonical_aa(std::uint8_t code) {
    return code < NUM_CANONICAL_AA;
  }

}

#endif // _PROTSEQCODE_HPP_
