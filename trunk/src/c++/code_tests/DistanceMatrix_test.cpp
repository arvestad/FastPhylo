//--------------------------------------------------
//                                        
// File: DistanceMatrix_test.cpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: DistanceMatrix_test.cpp,v 1.4 2006/09/04 11:54:18 isaac Exp $                                 
//
//--------------------------------------------------

#include "DistanceMatrix.hpp"
#include <fstream>
#include <iostream>
#include "file_utils.hpp"

using namespace std;

int
main(int argc, char **argv){

  std::ifstream in;
  open_read_stream("dm_test_file.txt",in);
  StrDblMatrix dm(in);
  in.close();
  
  cout << dm << endl;

  cout << "---------------" << endl;
  dm.setDistance(0,2,0);
  cout << dm << endl;

  cout << "---------------" << endl;
  StrDblMatrix dm2(dm);
  cout << dm2 << endl;


  cout << "--------------- SWAP " << endl;
  cout << dm2 << endl;
  PRINT_EXP(dm2.swapRowToLast(2));
  cout << dm2 << endl;

  PRINT_EXP(dm2.removeLastRow());
  cout << dm2 << endl;
  
  return 1;
}
