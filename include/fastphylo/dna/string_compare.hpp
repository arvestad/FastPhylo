//--------------------------------------------------
//
// File: string_compare.hpp
//
// Author: Isaac Elias
// e-mail: isaac@nada.kth.se
//
//--------------------------------------------------
#ifndef STRING_COMPARE_HPP
#define STRING_COMPARE_HPP

#include <array>
#include <string>
#include "fastphylo/dna/dna_pairwise_sequence_likelihood.hpp"
#include <iostream>

using divergence_matrix_t = std::array<std::array<float, 4>, 4>;

// Straigt forward hamming distance. Doesn't check for ambiguities or
// gaps.
int hamming_distance(const std::string &s1, const std::string &s2);

// Returns sufficient statistica for the Tamura Nei model
// i.e. the number of purine transitions, pyrimidine transitions, transversions, and the number of deletions
TN_string_distance TN_string_compare(const std::string &s1, const std::string &s2);

// returns the number of deleted characters.
int complete_dna_string_compare(
    divergence_matrix_t &divergence_matrix, // a 4x4 matrix will be filled in with the frequences
    const std::string &s1, const std::string &s2);

TN_string_distance divergence_matrix_2_TN_distance(const divergence_matrix_t &divergence_matrix, int deleted);

void print_divergence_matrix(std::ostream &out, const divergence_matrix_t &div_matrix);

#endif // STRING_COMPARE_HPP
