#pragma once

#include <vector>
#include "fastphylo/protein/ExpectedDistance.hpp" // for type_prior
#include "fastphylo/protein/ModelMatrix.hpp"      // for model_type
#include "fastphylo/protein/Matrix.hpp"           // for MatrixExpm
#include "fastphylo/core/DistanceMatrix.hpp"

// Forward Declarations
class Sequence;

using SeqVec = std::vector<Sequence>;

//! Struct that contains information needed for distance calculations
struct prot_sequence_translation_model
{
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
};

//! Calculate distances without standard deviation
void calculate_distances(const SeqVec &sv, StrDblMatrix &dm, prot_sequence_translation_model t_model);

//! Calculate distances with standard deviation
void calculate_distances(const SeqVec &sv, StrDblMatrix &dm, prot_sequence_translation_model t_model,
                         StrDblMatrix &sdm);

//! Calculates the expected distance without standard deviation
void calculate_ed_dists(const SeqVec &sv, StrDblMatrix &dm, const prot_sequence_translation_model &tm);

//! Calculates the expected distance with standard deviation
void calculate_ed_dists_with_sd(const SeqVec &sv, StrDblMatrix &dm, StrDblMatrix &sdm,
                                const prot_sequence_translation_model &tm, bool sd);

//! Builds the one-time decomposition calculate_ml_dists() needs for a
//! given named model - the expensive step (MatrixExpm's symmetrized
//! eigendecomposition, Matrix.hpp), so callers that invoke
//! calculate_ml_dists() more than once for the same model (bootstrap
//! replicates in fastprot/main.cpp; an a-priori-unknown number of
//! calls from a long-lived caller such as fastphylo-py bindings) must
//! call this ONCE and reuse the result, not call it again per
//! distance-estimation call - calculate_ml_dists() itself has no way
//! to build one from a model_type, by design, so it can't be redone
//! by accident inside a loop that only has access to that function.
MatrixExpm build_ml_decomposition(model_type mt);

//! Calculates the maximum likelihood distance for every pair in sv,
//! using an already-built decomposition - see build_ml_decomposition()
//! above for why this doesn't take a model_type and build one itself.
void calculate_ml_dists(const SeqVec &sv, StrDblMatrix &dm, const MatrixExpm &Qdecomp);

//! Calculates the identity distance
void calculate_id_dists(const SeqVec &sv, StrDblMatrix &dm);

//! Calculates the Jukes-Cantor corrected distance
void calculate_jc_dists(const SeqVec &sv, StrDblMatrix &dm);

//! Calculates the Kimura corrected distance
void calculate_kimura_dists(const SeqVec &sv, StrDblMatrix &dm);

//! Calculates the Storm-Sonnhammer corrected distance
void calculate_stormsonnhammer_dists(const SeqVec &sv, StrDblMatrix &dm);
