#pragma once

#include <vector>

// Forward Declarations
class Matrix;

using DblVec = std::vector<double>;

//! Enum that specifies the wanted model. The 10 models jtt_dcmut..pmb
//! were added 2026-08-11 to match fastphylo-py's own model catalog
//! (github.com/arvestad/fastphylo-py's src/fastphylo/protein.py) - see
//! ModelMatrix.cpp's table-driven get_model_matrix()/get_model_vec()
//! for how they're stored (differently from the original 5 - R +
//! equilibrium frequencies, not a literal precomputed Q, see that
//! file's comment for why). Appended, not alphabetized, so existing
//! enum integer values never shift.
enum model_type
{
    id,
    jc,
    jck,
    jcss,
    wag,
    day,
    jtt,
    mvr,
    lg,
    jtt_dcmut,
    vt,
    hivb,
    hivw,
    cprev,
    blosum62,
    dcmut,
    mtrev,
    rtrev,
    pmb
};

//! Gets the rate matrix for the specified model
Matrix get_model_matrix(model_type model);
//! Gets the equilibrium distribution for the specified model
DblVec get_model_vec(model_type model);
//! The mean substitution rate -sum(eq[i] * Q(i,i)) implied by Q at its
//! own, as-published scale, given its equilibrium distribution eq.
//! Rate matrices from different sources aren't necessarily published
//! at a consistent scale - e.g. WAG/JTT/Dayhoff/MVR are all ~100
//! (the historical "PAM" convention - expected substitutions per 100
//! residues), but LG is already published normalized to 1 (found via
//! cross-checking against IQ-TREE2's source, 2026-08-11: computed
//! directly from IQ-TREE2's own raw model/modelprotein.cpp data,
//! WAG/JTT/Dayhoff came out ~95-101, LG came out 1.000001). Dividing Q
//! by this value (see calculate_ml_dists(), ProtDistCalc.cpp) makes
//! `t` mean the same thing - expected substitutions/site - for every
//! model uniformly, regardless of which convention it was originally
//! published in, rather than assuming one fixed factor works for all
//! of them (it doesn't, for LG).
double mean_substitution_rate(const Matrix &Q, const DblVec &eq);
//! The JTT rate matrix
Matrix get_jtt();
//! The JTT equilibrium distribution
DblVec get_jtt_eq();
//! The Müller-Vingron rate matrix
Matrix get_mvr();
//! The Müller-Vingron equilibrium distribution
DblVec get_mvr_eq();
//! The Dayhoff rate matrix
Matrix get_day();
//! The Dayhoff equilibrium distribution
DblVec get_day_eq();
//! The WAG rate matrix
Matrix get_wag();
//! The WAG equilibrium distribution
DblVec get_wag_eq();
//! The LG rate matrix
Matrix get_lg();
//! The LG equilibrium distribution
DblVec get_lg_eq();
