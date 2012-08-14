#include "ProtSeqUtils.hpp"
#include <string>
#include <set>
#include "Matrix.hpp"
#include "../../Sequence.hpp"

/*
 * Function that translates an amino acid to a specific index
 * If the character given is not in the amino acid alphabet specified
 * below, 100 is return. Therefore the return value needs to be checked.
 * Amino acid order is: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V
 * @param c The aminoacid
 * @return The index
 */
  std::size_t getAAInd(char c){
    switch (c) {
      case 'A': return 0; 
      case 'R': return 1; 
      case 'N': return 2; 
      case 'D': return 3; 
      case 'C': return 4; 
      case 'Q': return 5; 
      case 'E': return 6; 
      case 'G': return 7; 
      case 'H': return 8; 
      case 'I': return 9; 
      case 'L': return 10; 
      case 'K': return 11; 
      case 'M': return 12; 
      case 'F': return 13; 
      case 'P': return 14; 
      case 'S': return 15; 
      case 'T': return 16; 
      case 'W': return 17; 
      case 'Y': return 18; 
      case 'V': return 19; 
      default: return 100;
    }
  }

  /* 
   * Saves the positions for all indels in the given sequences and
   * then erases all of them from every sequence.
   * @param sv A vector with sequences
   */
  void remove_gaps(std::vector<Sequence> &sv){
    std::set<int> positions;

    std::vector<Sequence>::iterator it;

    // Find all the gaps
    for (it = sv.begin(); it != sv.end(); it++){
      std::size_t found;

      found = it->seq.find_first_of('-');
      while (found != std::string::npos){
        positions.insert(found);
        found = it->seq.find_first_of('-', found+1);
      }
    }

    // Remove all the gaps
    for (it = sv.begin(); it != sv.end(); it++){
      std::set<int>::reverse_iterator set_rit;
      for (set_rit = positions.rbegin(); set_rit != positions.rend(); set_rit++){
        it->seq.erase(*set_rit, 1);
      }
    }
  }
/*
 * Counts all replacements from one amino acid to another in two sequences
 * @param s1 The first sequence
 * @param s2 The second sequence
 * @return A replacement count matrix
 */
Matrix count_replacements(const Sequence &s1, const Sequence &s2){
  Matrix temp(20,20);
  std::string::const_iterator it1, it2;
  for (it1 = s1.seq.begin(), it2 = s2.seq.begin(); it1 != s1.seq.end(); it1++, it2++) {
    int c1 = getAAInd(toupper(*it1));
    int c2 = getAAInd(toupper(*it2));
    if (c1 != 100 && c2 != 100){
      temp(c1, c2)++;
    }
  }
  return temp;
}

  /*
   * Counts the percentage identity, the number of amino acids that are 
   * the same between two sequences, divided with the length of the sequences
   * @param s1 Sequence 1
   * @param s2 Sequence 2
   * @return The percentage identity
   */
  double count_id_dist(const Sequence &s1, const Sequence &s2){
    double id = 0;
    std::string::const_iterator it1, it2;

    for (it1 = s1.seq.begin(), it2 = s2.seq.begin(); it1 != s1.seq.end(); it1++, it2++) {
      if (toupper(*it1) == toupper(*it2))
        id++;
    }
    return id/s1.seq.size();
  }

/*
 * Code adapted from Sequences2DistanceMatrix.cpp - bootstrapSequences()
 * @param seq Original vector with sequences
 * @param bseq Vector with new bootstrapped sequences
 */
void bootstrap_sequences(const std::vector<Sequence> &seqs, std::vector<Sequence> &bseqs){
  //ensure capacity in bootsequences
  bseqs.resize(seqs.size());

  const size_t seqlen = seqs[0].seq.length();
  for ( size_t i = 0 ; i < seqs.size() ; i++ )
    bseqs[i].seq.reserve(seqlen);

  // Do the bootstrapping
  size_t pos=0;
  size_t seq;  
  std::vector<int> samplePositions(seqlen);
  const size_t BUFFSIZE = (16383>seqlen ? seqlen : 16383); //2^14=16384
  char buff[BUFFSIZE+1];
  buff[BUFFSIZE]='\0';
  size_t const stride = 32;

  if( stride < seqlen){
    for( pos=0; pos<(seqlen-stride); pos+= stride)
      for( size_t i=0; i<stride; i++)
        samplePositions[pos+i] = (int) (seqlen*1.0*rand()/(RAND_MAX+1.0));
    for (; pos<seqlen; pos++)
      samplePositions[pos] = (int) (seqlen*1.0*rand()/(RAND_MAX+1.0));


    for( seq=0; seq<seqs.size(); seq++){
      const std::string & s = seqs[seq].seq;
      //-----------
      pos = 0;
      for( ; pos<seqlen-BUFFSIZE; pos+=BUFFSIZE){
        for( size_t i=0; i<BUFFSIZE; i++){
          buff[i] = s[samplePositions[pos+i]]; 
        }
        bseqs[seq].seq.append(buff);
      }
      size_t i;
      for ( i=0; pos<seqlen; i++,pos++)
        buff[i] = s[samplePositions[pos]]; 
      buff[i] = '\0';
      //----------
      bseqs[seq].seq.append(buff);
    }
  }
  else{//seqlen<stride
    for (pos=0; pos<seqlen; pos++)
      samplePositions[pos] = (int) (seqlen*1.0*rand()/(RAND_MAX+1.0));    

    for( seq=0; seq<seqs.size(); seq++){
      const std::string & s = seqs[seq].seq;
      for (pos=0; pos<seqlen; pos++)
        buff[pos] = s[samplePositions[pos]]; 
      buff[pos]='\0';
      bseqs[seq].seq.append(buff);
    } 
  }

}

