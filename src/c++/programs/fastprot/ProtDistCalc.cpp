#include "ProtDistCalc.hpp"
#include <cctype> // for toupper()
#include <cmath>  // for log() and pow()
#include <string>
#include "MaximumLikelihood.hpp"
#include "ProtSeqUtils.hpp"
#include "ProtSeqCode.hpp"
#include "ProtSeqCompare.hpp"
#include "Matrix.hpp"
#include "../../Sequence.hpp"

// ID, JC, JCK (Kimura), JCSS, ED and ML all reduce to a function of either
// the identity fraction or the replacement tally between two sequences.
// As of speed2026a these are always computed via the byte-encoded SIMD/
// encoded-array primitives (ProtSeqCode.hpp/ProtSeqCompare.hpp) - see
// phase0_audit.md/phase1_design.md/benchmarks/RESULTS.md for the
// profiling, design and benchmark work behind this, and plan.md for the
// original brief. The old scalar count_id_dist()/count_replacements()
// (toupper()+compare, or toupper()+getAAInd()+bounds-checked Matrix
// increments) have been validated byte-identical against their
// replacements and removed from this file.
//
// Note on scope: profiling ED (calculate_ed_dists_with_sd) and ML
// (calculate_ml_dists) found count_replacements() itself was only ~1.3%
// of ED's wall time (dominated by posterior_probability()'s per-pair
// elem_mult()+sum() allocation churn in ExpectedDistance.cpp) and
// didn't register at all in ML's profile (dominated by Matrix::expm()'s
// LAPACK eigendecomposition). Wiring in the faster primitive here is a
// real, safe win on the primitive itself (see benchmarks/RESULTS.md),
// not a fix for either of those larger, separate bottlenecks - neither
// is in scope for this change.
namespace {
  // Encodes every sequence in sv once (ProtSeqCode.hpp), so the per-pair
  // cost in the loops below doesn't include re-encoding.
  std::vector< std::vector<std::uint8_t> > encode_all(const SeqVec &sv) {
    std::vector< std::vector<std::uint8_t> > encoded(sv.size());
    for (std::size_t i = 0; i < sv.size(); i++)
      ProtSeqCode::encode_sequence(sv[i].seq, encoded[i]);
    return encoded;
  }

  // Converts ProtSeqCode::count_replacement_tally()'s flat
  // NUM_CANONICAL_AA x NUM_CANONICAL_AA tally into the Matrix type
  // ExpectedDistance.cpp/MaximumLikelihood.cpp expect, matching
  // count_replacements()'s existing (row,col) = (from,to) indexing.
  Matrix tally_to_matrix(const std::vector<std::size_t> &tally) {
    Matrix m(ProtSeqCode::NUM_CANONICAL_AA, ProtSeqCode::NUM_CANONICAL_AA);
    for (std::size_t a = 0; a < ProtSeqCode::NUM_CANONICAL_AA; a++)
      for (std::size_t b = 0; b < ProtSeqCode::NUM_CANONICAL_AA; b++)
        m(a, b) = static_cast<double>(tally[a * ProtSeqCode::NUM_CANONICAL_AA + b]);
    return m;
  }

  Matrix replacement_tally(const std::vector<std::uint8_t> &e1, const std::vector<std::uint8_t> &e2) {
    return tally_to_matrix(ProtSeqCode::count_replacement_tally(e1.data(), e1.size(), e2.data(), e2.size()));
  }
}

  /*
   * Method to call to calculate distances without standard deviation 
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   */
  void calculate_distances(const SeqVec &sv, StrDblMatrix &dm, prot_sequence_translation_model t_model){
    switch (t_model.model){
      case id: 
        calculate_id_dists(sv, dm); 
        break;
      case jc: 
        calculate_jc_dists(sv, dm); 
        break;
      case jck: 
        calculate_kimura_dists(sv, dm); 
        break;
      case jcss: 
        calculate_stormsonnhammer_dists(sv, dm); 
        break;
      case wag:
      case day:
      case arve:
      case jtt:
      case mvr:
      case lg:
        if (t_model.ml)
          calculate_ml_dists(sv, dm, t_model.model); 
        else 
          calculate_ed_dists_with_sd(sv, dm, dm, t_model, false); 
        break;
    }
  }

  /*
   * Method to call to calculate distances with standard deviation 
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   */
  void calculate_distances(const SeqVec &sv, StrDblMatrix &dm, 
      prot_sequence_translation_model t_model, StrDblMatrix &sdm){
    
    switch (t_model.model){
      case id: 
        calculate_id_dists(sv, dm); 
        break;
      case jc: 
        calculate_jc_dists(sv, dm); 
        break;
      case jck: 
        calculate_kimura_dists(sv, dm); 
        break;
      case jcss: 
        calculate_stormsonnhammer_dists(sv, dm); 
        break;
      case wag:
      case day:
      case arve:
      case jtt:
      case lg:
      case mvr:
        if (t_model.ml)
          calculate_ml_dists(sv, dm, t_model.model); 
        else 
          calculate_ed_dists_with_sd(sv, dm, sdm, t_model, true); 
        break;
    default:
      throw std::logic_error("Model not implemented");
    }

  }

  /*
   * Calculates the expected distance
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   * @param tm Translation model that contains parameters for the calculations
   */
  void calculate_ed_dists(const SeqVec &sv, StrDblMatrix &dm, 
      const prot_sequence_translation_model &tm){

    Matrix Q(get_model_matrix(tm.model));
    DblVec eq = get_model_vec(tm.model);
    initialize_ed(tm.tp, tm.step_size, Q, eq);
    dm.resize(sv.size());
    std::vector< std::vector<std::uint8_t> > encoded = encode_all(sv);

    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        Matrix N = replacement_tally(encoded[i], encoded[j]);
        double distance = 0.01 * calculate_ed(N);
        dm.setDistance(i, j, distance);
      }
    }

  }

  /*
   * Calculates the expected distance
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   * @param sdm Distance matrix where the standard deviation is saved, NULL if 
   * the standard deviation shouldn't be calculated
   * @param tm Translation model that contains parameters for the calculations
   */
  void calculate_ed_dists_with_sd(const SeqVec &sv, StrDblMatrix &dm, StrDblMatrix &sdm, 
      const prot_sequence_translation_model &tm, bool sd){

    Matrix Q = get_model_matrix(tm.model);
    DblVec eq = get_model_vec(tm.model);
    initialize_ed(tm.tp, tm.step_size, Q, eq);
    dm.resize(sv.size());
    if (sd)
      sdm.resize(sv.size());
    std::vector< std::vector<std::uint8_t> > encoded = encode_all(sv);

    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        Matrix N = replacement_tally(encoded[i], encoded[j]);
        double distance;
        if (sd)
          distance = 0.01 * calculate_ed_with_sd(N);
        else 
          distance = 0.01 * calculate_ed(N);
        dm.setDistance(i, j, distance); 
        
        if (sd)
          sdm.setDistance(i, j, get_standard_deviation());
      }
    }

  }
  /*
   * Calculates the maximum likelihood distance.
   *
   * Notice the rescaling with a factor 100. The Q matrices have a funny scaling due a bad decision years ago!
   *
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   * @param mt Specifies which model to use for the calculations
   */
  void calculate_ml_dists(const SeqVec &sv, StrDblMatrix &dm, model_type mt){
    Matrix Q = get_model_matrix(mt);
    // Q is the same rate matrix for every pair below; decomposing it is
    // the expensive part of evaluating exp(Q*t) (LAPACK eigendecomposition
    // + matrix inversion), so do it once here instead of once per Newton-
    // Raphson iteration per pair inside likelihood_calc()/
    // likelihood_deriv() - see Matrix.hpp's MatrixExpm and
    // phase0_audit.md's "ML speedup round" for the profiling behind this.
    MatrixExpm Qdecomp(Q);
    dm.resize(sv.size());
    std::vector< std::vector<std::uint8_t> > encoded = encode_all(sv);

    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        Matrix N = replacement_tally(encoded[i], encoded[j]);
        double distance = 0.01 * likelihood_calc(N, Q, Qdecomp);
        dm.setDistance(i, j, distance);
      }
    }
  }

  /*
   * Computes the identity based distance
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   */
  void calculate_id_dists(const SeqVec &sv, StrDblMatrix &dm){
    dm.resize(sv.size());
    std::vector< std::vector<std::uint8_t> > encoded = encode_all(sv);
    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        double distance = ProtSeqCode::count_id_fraction(encoded[i].data(), encoded[i].size(),
            encoded[j].data(), encoded[j].size());
        dm.setDistance(i, j, distance);
      }
    }
  }

  /*
   * Calculates the Jukes-Cantor corrected identity distance
   * d = -(19/20)*log(1-20/19)*(1-count_id_dist(s1, s2))
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   */
  void calculate_jc_dists(const SeqVec &sv, StrDblMatrix &dm){
    dm.resize(sv.size());
    std::vector< std::vector<std::uint8_t> > encoded = encode_all(sv);
    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        double diff = 1 - ProtSeqCode::count_id_fraction(encoded[i].data(), encoded[i].size(),
            encoded[j].data(), encoded[j].size());
        double distance = -(19/20.0)*log(1-(20.0/19) * diff);
        dm.setDistance(i, j, distance);
      }
    }

  }
  /*
   * Calculates the Kimura corrected distance between two sequences
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   */
  void calculate_kimura_dists(const SeqVec &sv, StrDblMatrix &dm){
    dm.resize(sv.size());
    std::vector< std::vector<std::uint8_t> > encoded = encode_all(sv);
    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        double diff = 1 - ProtSeqCode::count_id_fraction(encoded[i].data(), encoded[i].size(),
            encoded[j].data(), encoded[j].size());
        double adj_distance = diff + 0.2*diff*diff;
        if (adj_distance > 0.854)
          adj_distance = 0.854;
        double distance = - log(1-adj_distance);
        dm.setDistance(i, j, distance);
      }
    }
  }

  /*
   * Calculates Storm-Sonnhammer distance
   * @param sv Vector with protein sequences
   * @param dm Distance matrix where the results is saved
   */
void calculate_stormsonnhammer_dists(const SeqVec &sv, StrDblMatrix &dm){
  dm.resize(sv.size());
  std::vector< std::vector<std::uint8_t> > encoded = encode_all(sv);
  for (int i=0; i<sv.size(); i++) {
    for (int j=i+1; j<sv.size(); j++){
      double diff = 1 - ProtSeqCode::count_id_fraction(encoded[i].data(), encoded[i].size(),
          encoded[j].data(), encoded[j].size());
      double adj_distance;
      if (diff > 0.916)
        adj_distance = 1000.0;
      else {
        adj_distance = - log(1 - 0.95844*diff
            - 0.69957 * pow(diff, 2)
            + 2.4955 * pow(diff, 3)
            - 4.6353 * pow(diff, 4)
            + 2.8076 * pow(diff, 5));
      }
      dm.setDistance(i, j, adj_distance); 
    }
  }
}
