#include "ProtSeqUtils.hpp"
#include <string>
#include <set>
#include "../../Sequence.hpp"

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
  // getAAInd()/count_id_dist()/count_replacements() used to live here;
  // replaced by ProtSeqCode::count_id_fraction()/count_replacement_tally()
  // (ProtSeqCode.hpp/ProtSeqCompare.hpp) as of speed2026a Phase 6 (id) and
  // the count_replacements wiring round (replacements) - see
  // phase0_audit.md/phase1_design.md.

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

