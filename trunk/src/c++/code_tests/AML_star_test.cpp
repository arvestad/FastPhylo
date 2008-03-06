//--------------------------------------------------
//                                        
// File: AML_star_test.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: AML_star_test.cpp,v 1.3 2006/01/31 16:31:08 isaac Exp $                                 
//
//--------------------------------------------------

#include "AML_star.hpp"

#include <string>
#include <iostream>
#include "log_utils.hpp"

using namespace std;

int
main ( int argc, char ** argv ){


  string a = "aaaaaaaaac";
  string b = "aaaaaacccc";
  string c = "aaaaaaaagg";
  string center;

  PRINT(a);
  PRINT(b);
  PRINT(c);
  PRINT(find_AML_star(a,b,c,center, P_DISTANCE));
  PRINT(center);

  cout << "-------------" << endl;
  a = "aaatttggg";
  b = "aaacccaaa";
  c = "tttcccggg";

  PRINT(a);
  PRINT(b);
  PRINT(c);
  PRINT(find_AML_star(a,b,c,center, P_DISTANCE));
  PRINT(center);

  
  cout << "-------------" << endl;
  a = "aaatttaaaaaaaaaaaaaaa";
  b = "aaaaaatttaaaaaaaaaaaa";
  c = "tttaaaaaaaaaaaaaaaaaa";

  PRINT(a);
  PRINT(b);
  PRINT(c);
  PRINT(find_AML_star(a,b,c,center, P_DISTANCE));
  PRINT(center);

  
  cout << "-------------" << endl;
  cout << " Parsimony solution " << endl;
  a = "aaatttaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
  b = "aaaaaatttaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
  c = "tttaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";

  PRINT(a);
  PRINT(b);
  PRINT(c);
  PRINT(find_AML_star(a,b,c,center, P_DISTANCE));
  PRINT(center);

  
  cout << "-------------" << endl;
  cout << " Non parsimonious" << endl;
  a = "aattttttttttttttttttaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
  b = "ataaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
  c = "taaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";

  PRINT(a);
  PRINT(b);
  PRINT(c);
  PRINT(find_AML_star(a,b,c,center, P_DISTANCE));
  PRINT(center);
}
