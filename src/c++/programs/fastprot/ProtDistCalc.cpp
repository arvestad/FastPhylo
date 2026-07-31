#include "ProtDistCalc.hpp"
#include <cctype> // for toupper()
#include <cmath>  // for log() and pow()
#include <cstdlib> // for getenv()
#include <string>
#include "MaximumLikelihood.hpp"
#include "ProtSeqUtils.hpp"
#include "ProtSeqCode.hpp"
#include "ProtSeqCompare.hpp"
#include "Matrix.hpp"
#include "../../Sequence.hpp"

namespace {
  // Switch between the old (toupper()+compare, ProtSeqUtils.cpp) and new
  // (byte-encoded SIMD, ProtSeqCode.hpp/ProtSeqCompare.hpp) count_id_dist()
  // primitive, for the correctness/benchmark harness (Phase 4/5) to compare
  // against each other. Old path is the default until that harness signs
  // off (design decision 5 in plan.md). Only affects the four distance
  // models built on count_id_dist() (ID, JC, JCK/Kimura, JCSS); ED and ML
  // (built on count_replacements()) are untouched by this switch - not
  // wired this round, see phase0_audit.md.
  bool use_new_prot_primitives() {
    static const bool use = (std::getenv("FASTPROT_NEW_PROT_CODE") != NULL);
    return use;
  }

  // Encodes every sequence in sv once (ProtSeqCode.hpp), so the new
  // primitive's per-pair cost doesn't include re-encoding - matches how
  // Phase 5's benchmark needs to use it for a fair old-vs-new comparison.
  std::vector< std::vector<std::uint8_t> > encode_all(const SeqVec &sv) {
    std::vector< std::vector<std::uint8_t> > encoded(sv.size());
    for (std::size_t i = 0; i < sv.size(); i++)
      ProtSeqCode::encode_sequence(sv[i].seq, encoded[i]);
    return encoded;
  }

  // Identity fraction via whichever primitive is selected. encoded is NULL
  // when the old path is selected (encode_all() wasn't run).
  double id_fraction(const SeqVec &sv, const std::vector< std::vector<std::uint8_t> > *encoded,
      std::size_t i, std::size_t j) {
    if (encoded)
      return ProtSeqCode::count_id_fraction((*encoded)[i].data(), (*encoded)[i].size(),
          (*encoded)[j].data(), (*encoded)[j].size());
    return count_id_dist(sv[i], sv[j]);
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

    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        Matrix N = count_replacements(sv[i], sv[j]);
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

    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        Matrix N = count_replacements(sv[i], sv[j]);
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
    dm.resize(sv.size());

    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        Matrix N = count_replacements(sv[i], sv[j]);
        double distance = 0.01 * likelihood_calc(N, Q);
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
    bool use_new = use_new_prot_primitives();
    std::vector< std::vector<std::uint8_t> > encoded;
    if (use_new) encoded = encode_all(sv);
    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        double distance = id_fraction(sv, use_new ? &encoded : NULL, i, j);
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
    bool use_new = use_new_prot_primitives();
    std::vector< std::vector<std::uint8_t> > encoded;
    if (use_new) encoded = encode_all(sv);
    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        double diff = 1 - id_fraction(sv, use_new ? &encoded : NULL, i, j);
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
    bool use_new = use_new_prot_primitives();
    std::vector< std::vector<std::uint8_t> > encoded;
    if (use_new) encoded = encode_all(sv);
    for (int i=0; i<sv.size(); i++) {
      for (int j=i+1; j<sv.size(); j++){
        double diff = 1 - id_fraction(sv, use_new ? &encoded : NULL, i, j);
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
  bool use_new = use_new_prot_primitives();
  std::vector< std::vector<std::uint8_t> > encoded;
  if (use_new) encoded = encode_all(sv);
  for (int i=0; i<sv.size(); i++) {
    for (int j=i+1; j<sv.size(); j++){
      double diff = 1 - id_fraction(sv, use_new ? &encoded : NULL, i, j);
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
