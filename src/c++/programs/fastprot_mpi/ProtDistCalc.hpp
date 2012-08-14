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
  //! Speed
  int step_size;
  //! What kind of prior probability to be used, normal or flat
  type_prior tp;
} prot_sequence_translation_model;


  //! Calculate distances 
  void calculate_distances(const SeqVec &sv, StrDblMatrix &dm, 
      prot_sequence_translation_model t_model); 
  
  //! Calculates the expected distance 
  void calculate_ed_dists(const SeqVec &sv, StrDblMatrix &dm, 
      const prot_sequence_translation_model &tm);
 
  //! MPI version of the previous ED function
  void calculate_ed_dists_mpi(const SeqVec &sv, std::vector<double> &dv,
	const prot_sequence_translation_model &tm, int rank, int size);
 
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
