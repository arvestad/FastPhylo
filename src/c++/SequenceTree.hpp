//--------------------------------------------------
//                                        
// File: SequenceTree.hpp                             
//                             
// Author: Mehmood Alam Khan, Isaac Elias
// e-mail: malagori@kth.se, isaac@nada.kth.se
//                             
// cvs: $Id: SequenceTree.hpp,v 1.22 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------
#ifndef SEQUENCETREE_HPP
#define SEQUENCETREE_HPP

#include "Tree.hpp"
#include <string>
#include "InitAndPrintOn_utils.hpp"
#include "Sequence.hpp"
#include "DistanceMatrix.hpp"
#include "FloatDistanceMatrix.hpp"
#include "BitVector.hpp"

//---------------------------------------
// A REGULAR SEQUENCE TREE
// Each node has a struct of a sequence and a double.
//--------------------------------------

//access defines for the node data.
#define NAME(node) (node)->data.s.name
#define SEQ(node) (node)->data.s.seq
#define SEQUENCE(node) (node)->data.s
#define EDGE(node) (node)->data.dbl





class SequenceTree : public Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double> >
{
public:
  //the same constructers as in Tree
  SequenceTree() : Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double> >(){}
  
  SequenceTree(Sequence_double d) : Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double> >(d){}
  
  SequenceTree(const SequenceTree &t) : 
    Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double> >(t){}
  
  SequenceTree(const TreeNode<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double> > &n) :
    Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double> >(n){}
  
  SequenceTree(char *newickstr) : 
    Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double> >(newickstr){}
  SequenceTree(std::istream &in) :
    Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double> >(in){}

  template<class Data2, class DataInit2, class DataPrintOn2>
  SequenceTree(const Tree<Data2,DataInit2,DataPrintOn2> &t, Sequence_double defaultData) :
    Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double> >(t,defaultData){}
  
  virtual ~SequenceTree(){}

  

  //----------------------------------
  // ADDITION PRINT FUNCTIONS
  // The print functions in Tree are not shadowed.
  void verbosePrint(std::ostream &os);
  void printNodeData(std::ostream &os);
  //sets the names on nodes without names to be "n"+nodeId
  void setNodeNames();
  static void printSequencesPhylip(SequenceTree::NodeVector &nodes,std::ostream &os);
  void printSequencesPhylip(std::ostream &os);
  void printSequences(std::ostream &os);
  void printSequencesWithoutGaps(std::ostream &os);


  
  //------------------
  //removes all degree 2 nodes and makes an edge length a+b where
  //a and b are the lengths of the the two edges
  void shortcutDegree2Nodes();
  void shortcutFloatDegree2Nodes();

  //---------------
  //returns the sum of all edge lengths
  double sumOfEdgeLengths();
  float sumOfFloatEdgeLengths();


  //---------------
  //COMPUTE LIKELIHOOD
  //Computes the likelihood of the tree.  If the an edge length is
  //negative it is chosen to be the optimal edge length.
  // This function works for p-distance and JukesCantor
  // since the likelihood of these are the same
  // All internal nodes have to have strings assigned to them.
  double compute_loglikelihood();
  float computeFloat_loglikelihood();

  //--------------------------- 
  // MAKE CANONICAL
  //the leafs are added in order of their ids in name2id, returns false
  //if the name2id map doesn't contain all leafs
  void makeCanonical(const str2int_hashmap &name2id);
  //the makeCanonical in the super class
  using Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double> >::makeCanonical;

  //takes a tree and fills the name2id map
  //such that each leaf name gets an id in the range [0,numleafs)
  void createLeafNameToLeafIdMap(str2int_hashmap &name2id) const;
  

  //-------------------
  // ROBINSON-FOULDS
  // takes time O(n^2) (this can be done in linear time)
  //i.e RF(T1,T2) = |Splitts(T1)\Splitts(T2)| + |Splitts(T2)\Splitts(T1)|
  //The normalized meassure is thus RF(T1,T2)/(Splitts(T1)+Splitts(T2))
  static double
  computeRobinsonFoulds(SequenceTree &t1,SequenceTree &t2);
  
  static float
    computeFloatRobinsonFoulds(SequenceTree &t1,SequenceTree &t2);

  void
  computeSplittSet(std::vector<BitVector> &splitts, SequenceTree::NodeVector &nodeForEachSplitt,
		   str2int_hashmap &name2index);
  static void
  printSplitt(BitVector &splitt, str2int_hashmap &name2index, std::ostream& os=std::cout);

  

  //------------ 
  //Takes a tree with all internal sequences set and
  //computes the hamming distance for each edge.
  void computeEdgeLengths();

  //PARSIMONY 
  //computes the most parsimonious sequences given the
  //sequences at the leafs.
  size_t computeMostParsimoniousSequences();
  
  //-----------------------------
  //
  // A distance matrix with identifiers as tree nodes.
  // Used by NJ and Big_AML
  typedef DistanceMatrix<SequenceTree::Node *,double,
			 Data_init<SequenceTree::Node *>,Data_printOn<SequenceTree::Node *>,
			 Data_init<double>,Data_printOn<double> > NodeMatrix;
  
  void tree2distanceMatrix(StrDblMatrix &dm);
  

  typedef FloatDistanceMatrix<SequenceTree::Node *,float,
  			 Data_init<SequenceTree::Node *>,Data_printOn<SequenceTree::Node *>,
  			 Data_init<float>,Data_printOn<float> > NodeFloatMatrix;

  void tree2FloatdistanceMatrix(StrFloMatrix &fdm);
  //--------------------------------------
  //All edges <=bound are removed.
  //The number of contracted edges are returned.
  //Edges incident to leafs are not removed. 
  int
  contractEdgesShorterThan(double bound=0);
  int
   contractFloatEdgesShorterThan(float bound=0);
  //--------
  //takes an array of strings which should consist of name,seq pairs/
  //i.e. map{n1,s1,n2,s2,...nN,sN} 
  void
  mapSequencesOntoTree( char  **nameseqPairs, int numPairs);

  void
  mapSequencesOntoTree(std::vector<Sequence> &seqs);

  //Takes a phylip sequence file and maps the sequences onto the tree
  //nodes by looking up the name.
  //It is possible to do the mapping from several files. That is some names
  //may exist in one file and other names in another file.
  void
  mapSequencesOntoTree(std::istream &in);

};

typedef __gnu_cxx::hash_map<const SequenceTree , int, objhash, objeq> tree2int_map;

#endif // SEQUENCETREE_HPP









