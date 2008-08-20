#ifndef NEWICK_HPP
#define NEWICK_HPP
#include <string>

struct NewickDelimiters {
   std::string left_parenthesis;
   std::string right_parenthesis;
   std::string left_leaf;
   std::string right_leaf;
   std::string comma;
   std::string null_tree;
   std::string semi_colon;
   std::string sequence_double_left;
   std::string sequence_double_right;

  NewickDelimiters() : left_parenthesis("("), right_parenthesis(")"), comma(","), semi_colon(";"), left_leaf(""), right_leaf(""),
		       sequence_double_left(":"),  sequence_double_right(""), null_tree("NULLTREE") {}
};


extern NewickDelimiters newickDelimiters;
#endif // NEWICK_HPP
