//--------------------------------------------------
//                                        
// File: Leastsquaresfit.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: LeastSquaresFit.hpp,v 1.3 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef LEASTSQUARESFIT_HPP
#define LEASTSQUARESFIT_HPP

#include "SequenceTree.hpp"
#include "DistanceMatrix.hpp"
#include <string>

#include "log_utils.hpp"
#include <iostream>
#include <float.h>

//-------------------------------------------------------------
//Takes a distance matrix in which the identifiers are tree nodes.
//and computes the least squares fit of all the edge lengths.
//Returns the least squares error.




double
computeLeastSquaresEdgeLengths(const StrDblMatrix &orig_dm,  SequenceTree &resultTree);

// computes the L2 distance betweent he two matrices.
double
computeL2(const StrDblMatrix &A,  const StrDblMatrix &B);

#endif // LEASTSQUARESFIT_HPP










