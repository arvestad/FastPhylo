#ifndef _MODEL_MAT_HPP_
#define _MODEL_MAT_HPP_

#include <vector>

  // Forward Declarations
  class Matrix;

  typedef std::vector<double> DblVec;
  
  //! Enum that specifies the wanted model
  enum model_type {id, jc, jck, jcss, wag, day, arve, jtt, mvr};

  //! Gets the rate matrix for the specified model
  Matrix get_model_matrix(model_type model);
  //! Gets the equilibrium distribution for the specified model
  DblVec get_model_vec(model_type model);
  //! The Arvestad rate matrix
  Matrix get_arvestad();
  //! The Arvestad equilibrium distribution
  DblVec get_arve_eq();
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
#endif
