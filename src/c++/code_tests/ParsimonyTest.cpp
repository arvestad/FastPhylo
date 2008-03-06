#include "SequenceTree.hpp"



using namespace std;

int
main(int argc,
     char **argv){

  
  SequenceTree tree("( (a:2,b:3):4  , (e:5,f:6):8)");
  
  SequenceTree::NodeVector nodes;
  tree.addLeafs(nodes);
  

  SEQ(nodes[0]) = "aaa";
  SEQ(nodes[1]) = "aag";
  SEQ(nodes[2]) = "aca";
  SEQ(nodes[3]) = "aac";

  PRINT(tree);
  tree.printSequencesPhylip(cout);

  PRINT(tree.computeMostParsimoniousSequences());

  
  PRINT(tree);
  tree.printSequencesPhylip(cout);

  return 1;    
  
}
