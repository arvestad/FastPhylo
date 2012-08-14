#ifndef _PROTDISTCALC_HPP
#define _PROTDISTCALC_HPP

#include <vector>
#include "ExpectedDistance.hpp"        // for type_prior
#include "ModelMatrix.hpp"             // for model_type
#include "../../DistanceMatrix.hpp"

  // Forward Declarations
class Sequence;
class Matrix;

typedef std::vector<Sequence> SeqVec;

//! Struct that contains information needed for distance calculations
typedef struct {
  //! Specifies what model to use
  model_type model;
  //! If a maximum likelihood computation should be made
  bool ml;
  //! If the standard deviations should be calculated
  bool sd;
  //! Speed
  int step_size;
  //! What kind of prior probability to be used, normal or flat
  type_prior tp;
} prot_sequence_translation_model;


  //! Calculate distances without standard deviation
  void calculate_distances(const SeqVec &sv, StrDblMatrix &dm, 
      prot_sequence_translation_model t_model); 
  
  //! Calculate distances with standard deviation
  void calculate_distances(const SeqVec &sv, StrDblMatrix &dm, 
      prot_sequence_translation_model t_model, StrDblMatrix &sdm); 
  
  //! Calculates the expected distance without standard deviation
  void calculate_ed_dists(const SeqVec &sv, StrDblMatrix &dm, 
      const prot_sequence_translation_model &tm);
  
  //! Calculates the expected distance with standard deviation
  void calculate_ed_dists_with_sd(const SeqVec &sv, StrDblMatrix &dm, StrDblMatrix &sdm, 
      const prot_sequence_translation_model &tm, bool sd);
  
  //! Calculates the maximum likelihood distance
  void calculate_ml_dists(const SeqVec &sv, StrDblMatrix &dm, model_type mt);
  
  //! Calculates the identity distance
  void calculate_id_dists(const SeqVec &sv, StrDblMatrix &dm);
  
  //! Calculates the Jukes-Cantor corrected distance
  void calculate_jc_dists(const SeqVec &sv, StrDblMatrix &dm);
  
  //! Calculates the Kimura corrected distance
  void calculate_kimura_dists(const SeqVec &sv, StrDblMatrix &dm);
  
  //! Calculates the Storm-Sonnhammer corrected distance
  void calculate_stormsonnhammer_dists(const SeqVec &sv, StrDblMatrix &dm);

#endif
