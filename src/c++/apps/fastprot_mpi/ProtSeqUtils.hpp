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
  
  //! Calculates percentage identity
  double count_id_dist(const Sequence &s1, const Sequence &s2);
  
  
  //! Performs bootstrapping
  void bootstrap_sequences(const std::vector<Sequence> &seqs, std::vector<Sequence> &bseqs); 
#endif
