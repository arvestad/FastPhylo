
#include "LeastSquaresFit.hpp"
#include "log_utils.hpp"
#include <iostream>

using  namespace std;
int
main(int argc,
     char **argv){
  

  SequenceTree tree("((a:1,b:2):3,c:2,d:1);");
  
  StrDblMatrix dm(tree.getNumLeafs());
  tree.tree2distanceMatrix(dm);
  tree.drawTree(cout);
  cout << dm <<endl;
  
  SEPARATOR();
  dm.setDistance(0,1,2);
  cout << dm << endl;

  double fit = computeLeastSquaresEdgeLengths(dm,tree);
  PRINT(fit);
  ASSERT_EQ(fit,0);
  tree.drawTree(cout);
  StrDblMatrix dm2(tree.getNumLeafs());
  tree.tree2distanceMatrix(dm2);
  cout << dm2 << endl;
  return 1;
}





