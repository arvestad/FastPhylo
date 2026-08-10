#pragma once

#include <vector>

// Forward Declarations
class Matrix;

using DblVec = std::vector<double>;

//! Enum that specifies the wanted model
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
    lg
};

//! Gets the rate matrix for the specified model
Matrix get_model_matrix(model_type model);
//! Gets the equilibrium distribution for the specified model
DblVec get_model_vec(model_type model);
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
