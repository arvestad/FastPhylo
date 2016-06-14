//--------------------------------------------------
//                                        
// File: SequenceTree.cpp                              
//                             
// Author: Mehmood Alam Khan,Isaac Elias
// e-mail: malagori@kth.se, isaac@nada.kth.se
//                             
// cvs: $Id: SequenceTree.cpp,v 1.44 2006/12/08 11:09:13 isaac Exp $                                 
//
//--------------------------------------------------

#include "SequenceTree.hpp"

#include <vector>
#include <string>
#include <unordered_set>

#include "InitAndPrintOn_utils.hpp"
#include <iostream>
#include "string_compare.hpp"
#include "stl_utils.hpp"
#include "log_utils.hpp"
#include "nucleotide.hpp"
#include "BitVector.hpp"

using namespace std;

//=---------------------
void
SequenceTree::verbosePrint(std::ostream &os){
  os << "---- SequenceTree " << endl;
  os << (*this) << endl << endl;
  drawTree(os);
  os << endl;
  printNodeData(os);
  os <<"----" <<endl;
}


void
SequenceTree::printNodeData(std::ostream &os){
  
  SequenceTree::NodeVector vec;
  vec.reserve(getNumNodes());
  addNodesInInfixOrder(vec);
  for ( size_t i = 0 ; i < vec.size() ; i++ ){
    string str = vec[i]->data.s.toString();
    if ( str.size() > 0 )
      os << "[" << vec[i]->getNodeId() << "] " << str << endl;
  }
}

void
SequenceTree::setNodeNames(){
  SequenceTree::NodeVector vec;
  vec.reserve(getNumNodes());
  addNodesInInfixOrder(vec);
  for ( size_t i = 0 ; i < vec.size() ; i++ ){
    SequenceTree::Node * n = vec[i];
    if ( NAME(n).length() == 0 )
      NAME(n) = string("n")+n->getNodeId();
  } 
}
void
SequenceTree::printSequencesPhylip(SequenceTree::NodeVector &nodes,  std::ostream &os){
  os << nodes.size() << "\t " << SEQ(nodes[0]).length() << endl;
  for ( size_t i = 0 ; i < nodes.size() ; i++ ){
    SequenceTree::Node *n = nodes[i];
    os << n->data.s << endl;
  } 
}
void 
SequenceTree::printSequencesPhylip(std::ostream &os){
  SequenceTree::NodeVector nodes;
  addNodesInInfixOrder(nodes);
  printSequencesPhylip(nodes,os);
}
void
SequenceTree::printSequences( std::ostream &os){
  SequenceTree::NodeVector vec;
  vec.reserve(getNumNodes());
  addNodesInInfixOrder(vec);
  for ( size_t i = 0 ; i < vec.size() ; i++ ){
    SequenceTree::Node *n = vec[i];
    if ( NAME(n).length() == 0 )
      continue;
    os << n->data.s << endl;
  } 
}
void 
SequenceTree::printSequencesWithoutGaps(std::ostream &os){
  SequenceTree::NodeVector vec;
  vec.reserve(getNumNodes());
  addNodesInInfixOrder(vec);
  for ( size_t i = 0 ; i < vec.size() ; i++ ){
    SequenceTree::Node *n = vec[i];
    if ( NAME(n).length() == 0 )
      continue;
    os << n->data.s << endl;
  } 
}


double
SequenceTree::sumOfEdgeLengths(){
  SequenceTree::NodeVector nodes;

  addNodesInPrefixOrder(nodes);
  double sum =0;
  //skip the root
  for(size_t i=1;i<nodes.size();i++)
    sum+=EDGE(nodes[i]);

  return sum;
}

float
SequenceTree::sumOfFloatEdgeLengths(){
  SequenceTree::NodeVector nodes;

  addNodesInPrefixOrder(nodes);
  float sum =0;
  //skip the root
  for(size_t i=1;i<nodes.size();i++)
    sum+=EDGE(nodes[i]);

  return sum;
}
//--------------------
//SHORTCUT
void SequenceTree::shortcutDegree2Nodes(){
  SequenceTree::NodeVector vec;
  addNodesWithDegree(vec,2);

  for ( size_t i = 0 ; i < vec.size() ; i++ ){
    SequenceTree::Node *n = vec[i];
    if ( vec[i]->isRoot() ){
      SequenceTree::Node *child = n->getRightMostChild();
      double newedge = EDGE(child) + EDGE(child->getLeftSibling());
      shortcutNode(vec[i]);
      EDGE(child) = -1;
      EDGE(child->getRightMostChild()) = newedge;
    }
    else{
      double pedge = EDGE(n);
      EDGE(n->getRightMostChild()) += pedge;
      shortcutNode(vec[i]);
    }
  }
  if ( getRoot()->getDegree() == 1 ){
    reRootAtHighDegreeNode();
    SequenceTree::Node *root = getRoot();
    EDGE(root) = -1;
  }
}

void SequenceTree::shortcutFloatDegree2Nodes(){
  SequenceTree::NodeVector vec;
  addNodesWithDegree(vec,2);

  for ( size_t i = 0 ; i < vec.size() ; i++ ){
    SequenceTree::Node *n = vec[i];
    if ( vec[i]->isRoot() ){
      SequenceTree::Node *child = n->getRightMostChild();
      float newedge = EDGE(child) + EDGE(child->getLeftSibling());
      shortcutNode(vec[i]);
      EDGE(child) = -1;
      EDGE(child->getRightMostChild()) = newedge;
    }
    else{
      float pedge = EDGE(n);
      EDGE(n->getRightMostChild()) += pedge;
      shortcutNode(vec[i]);
    }
  }
  if ( getRoot()->getDegree() == 1 ){
    reRootAtHighDegreeNode();
    SequenceTree::Node *root = getRoot();
    EDGE(root) = -1;
  }
}


//---------------
//COMPUTE LIKELIHOOD

double
SequenceTree::compute_loglikelihood(){

  SequenceTree::NodeVector vec;
  addNodesInPrefixOrder(vec);

  double hamdist;
  
  //skip the root node
  double loglikelihood = 0;
  for ( size_t i = 1 ; i < vec.size() ; i++ ){
    SequenceTree::Node *n = vec[i];
    double slen = 1.0*SEQ(n).size();
    hamdist = hamming_distance(SEQ(n),SEQ(n->getParent()));
    //   PRINT(n->data.s.name); PRINT(n->getParent()->data.s.name);PRINT(hamdist);
    //cout <<"-----------"<< endl;
    // if ( n->data.flt >=0 ){
    //       //cout << "prob: " << n->data.flt << endl;
    //       if ( hamdist != 0 ){
    //         PRINT(hamdist);
    //         PRINT(n->data.flt);
    //         PRINT(hamdist*log(n->data.flt/3.0) + (slen - hamdist)*log(1 - n->data.flt));
    //         loglikelihood += hamdist*log(n->data.flt/3.0) + (slen - hamdist)*log(1 - n->data.flt);
    //         PRINT(loglikelihood);
    //       }
    //     }
    //     else {//if no edge length available take optimadl edge length
    double optprob = hamdist/slen;
    optprob = ( 0.499 < optprob ? 0.499 : optprob );//we don't allow for more than 0.5 probability of change
    //cout << "optprob: " << optprob << endl;
    if ( hamdist != 0 ){
      //PRINT(hamdist*log(optprob/3.0) + (slen - hamdist)*log(1 - optprob));
      loglikelihood += hamdist*log(optprob/3.0) + (slen - hamdist)*log(1 - optprob);
      EDGE(n) = optprob;
    }
    else
      EDGE(n) = 0;
    //}
    //cout << "newlog " << loglikelihood << endl;
  }
  //cout << "------"<< endl;
  return loglikelihood;
}
float
SequenceTree::computeFloat_loglikelihood(){

  SequenceTree::NodeVector vec;
  addNodesInPrefixOrder(vec);

  float hamdist;

  //skip the root node
  float loglikelihood = 0;
  for ( size_t i = 1 ; i < vec.size() ; i++ ){
    SequenceTree::Node *n = vec[i];
    float slen = 1.0*SEQ(n).size();
    hamdist = hamming_distance(SEQ(n),SEQ(n->getParent()));
    //   PRINT(n->data.s.name); PRINT(n->getParent()->data.s.name);PRINT(hamdist);
    //cout <<"-----------"<< endl;
    // if ( n->data.flt >=0 ){
    //       //cout << "prob: " << n->data.flt << endl;
    //       if ( hamdist != 0 ){
    //         PRINT(hamdist);
    //         PRINT(n->data.flt);
    //         PRINT(hamdist*log(n->data.flt/3.0) + (slen - hamdist)*log(1 - n->data.flt));
    //         loglikelihood += hamdist*log(n->data.flt/3.0) + (slen - hamdist)*log(1 - n->data.flt);
    //         PRINT(loglikelihood);
    //       }
    //     }
    //     else {//if no edge length available take optimadl edge length
    float optprob = hamdist/slen;
    optprob = ( 0.499 < optprob ? 0.499 : optprob );//we don't allow for more than 0.5 probability of change
    //cout << "optprob: " << optprob << endl;
    if ( hamdist != 0 ){
      //PRINT(hamdist*log(optprob/3.0) + (slen - hamdist)*log(1 - optprob));
      loglikelihood += hamdist*log(optprob/3.0) + (slen - hamdist)*log(1 - optprob);
      EDGE(n) = optprob;
    }
    else
      EDGE(n) = 0;
    //}
    //cout << "newlog " << loglikelihood << endl;
  }
  //cout << "------"<< endl;
  return loglikelihood;
}

void
SequenceTree::computeEdgeLengths(){
  SequenceTree::NodeVector nodes;
  addNodesInPrefixOrder(nodes);
  for(size_t i=1 ; i<nodes.size() ; i++){
    EDGE(nodes[i]) = hamming_distance(SEQ(nodes[i]),SEQ(nodes[i]->getParent()));
  }
}


int
SequenceTree::contractEdgesShorterThan(double bound){
  SequenceTree::NodeVector nodes;
  addNodesInPrefixOrder(nodes);
  int numC=0;
  for(size_t i=1 ; i<nodes.size() ; i++){
    if(nodes[i]->isLeaf()) continue;
    if(EDGE(nodes[i])<=bound){
      shortcutNode(nodes[i]);
      numC++;
    }
  }
  return numC;
}
int
SequenceTree::contractFloatEdgesShorterThan(float bound){
  SequenceTree::NodeVector nodes;
  addNodesInPrefixOrder(nodes);
  int numC=0;
  for(size_t i=1 ; i<nodes.size() ; i++){
    if(nodes[i]->isLeaf()) continue;
    if(EDGE(nodes[i])<=bound){
      shortcutNode(nodes[i]);
      numC++;
    }
  }
  return numC;
}
//---------------------
typedef std::unordered_map<const std::string, SequenceTree::Node *, hashstr, eqstr> str2node_map;

void
SequenceTree::mapSequencesOntoTree(char  **nameseqPairs, int numPairs){

  size_t numnodes = getNumNodes();
  SequenceTree::NodeVector nodes;
  nodes.reserve(numnodes);
  addNodesInPrefixOrder(nodes);

    
  //build hash map {name,node}
  str2node_map str2node((int)(numnodes*1.5));
  for ( size_t i = 0 ; i < numnodes ; i++ ){
    if ( NAME(nodes[i]).size() > 0 ){
      str2node[NAME(nodes[i])] = nodes[i];
    }
  } 

  //go through the nameseq pairs and look up node
  for ( int i = 0 ; i < 2*numPairs ; i+=2 ){
    str2node_map::iterator iter = str2node.find(string(nameseqPairs[i]));
    if ( iter != str2node.end() ){
      SEQ((*iter).second).clear();
      SEQ((*iter).second).append(nameseqPairs[i+1]);
    }
    else{
      //      USER_WARNING("Unknown name in tree: " << nameseqPairs[i]);
    }
  }
}

void
SequenceTree::mapSequencesOntoTree( std::vector<Sequence> &seqs){
  size_t numnodes = getNumNodes();
  SequenceTree::NodeVector nodes;
  nodes.reserve(numnodes);
  addNodesInPrefixOrder(nodes);

    
  //build hash map {name,node}
  str2node_map str2node((int)(numnodes*1.5));
  for ( size_t i = 0 ; i < numnodes ; i++ ){
    if ( NAME(nodes[i]).size() > 0 ){
      str2node[NAME(nodes[i])] = nodes[i];
    }
  } 

  //go through the nameseq pairs and look up node
  for ( size_t i=0 ; i<seqs.size() ; i++ ){
    str2node_map::iterator iter = str2node.find(seqs[i].name);
    if ( iter != str2node.end() ){
      SEQ((*iter).second).clear();
      SEQ((*iter).second).append(seqs[i].seq);
    }
    else{
      //      USER_WARNING("Unknown name in tree: " << nameseqPairs[i]);
    }
  }

}

//--------------------
void
SequenceTree::mapSequencesOntoTree(std::istream &fin){

  int numSequences;
  unsigned int seqlen;
  char tmp[100];
  fin.getline(tmp,100);
  sscanf(tmp,"%d %d",&numSequences,&seqlen);
  

  //build hash map {name,node} and reserve mem for the sequences
  str2node_map str2node((int)(getNumNodes()*1.5));
  SequenceTree::NodeVector nodes;
  nodes.reserve(getNumNodes());
  addNodesInPrefixOrder(nodes);
  for ( size_t i = 0 ; i < nodes.size() ; i++ ){
    SEQ(nodes[i]).clear();
    SEQ(nodes[i]).reserve(seqlen+10);
    if ( NAME(nodes[i]).size() > 0 ){
      str2node[NAME(nodes[i])] = nodes[i];
    }
  }

  //names in the file are matched agains a node sequences If the name
  //doesn't exist in the tree then these sequences are read into a
  //garbage string.
  string garbage;
  string *actualNodeString = NULL;//used to check that the whole sequence has been read
  garbage.reserve(seqlen+10);
  std::vector<string> names(numSequences);
  std::vector<string *> sequences(numSequences,&garbage);

  
  //read the names and map the sequences onto the tree
  
  for ( int i = 0 ; i < numSequences ; i++ ){
    fin >> names[i];
    garbage.clear();
    
    //find the node 
    str2node_map::iterator iter = str2node.find(string(names[i]));
    if ( iter != str2node.end() ){
      sequences[i] = & SEQ(((*iter).second));
      actualNodeString = sequences[i];
    }

    string *currseq = sequences[i];
    
    while (1){
      char c = fin.get();
      nucleotide n = char2nucleotide(c);
      if ( DNA_NOT_ALLOWED == n ){
        if ( !isspace(c) ){
          USER_ERROR("Bad character \'" << c << "\'");
        }
        else if ( c != '\n' )
          continue;
        //if '\n'
        break;
      }

      currseq->append(1,nucleotide2char(n));
    }
    
  }
  
  //read remaining sequences
  
  //The sequences aren't neccesarily on one line but my be spread out interleaving
  //over several lines. Therefore we read until seqlen chars have been read.
  while ( actualNodeString->length() < seqlen ){
    for ( int i = 0 ; i < numSequences ; i++ ){
      char c = fin.peek();
      if ( !isspace(c) ){
        //skip first 10 chars
        for ( int j = 10 ; j != 0 ; j-- )
          fin.get();
      }
      while ( c == '\n' )
        c = fin.get();

      garbage.clear();
      string *currseq = sequences[i];
      while (1){
        c = fin.get();
        nucleotide n = char2nucleotide(c);
        if ( DNA_NOT_ALLOWED == n ){
          if ( !isspace(c) ){
            USER_ERROR("Bad character \'" << c << "\'");
          }
          else if ( c != '\n' )//skip space
            continue;
          //if '\n'
          break;
        }

        currseq->append(1,nucleotide2char(n));
      } 
    }
  }


  // CHECK THAT ALL STRINGS HAVE THE SAME LENGTH
  for ( size_t i = 0 ; i < nodes.size() ; i++ )
    if ( SEQ(nodes[i]).length() != seqlen && SEQ(nodes[i]).length() != 0 ){
      USER_ERROR("Sequence not of correct length: " << nodes[i]->data.s.name ); 
    }

}




//----------------
//-------------------
// ROBINSON-FOULDS

typedef std::unordered_set<BitVector*, objhash_ptr, objeq_ptr> BitVectorPtr_set;

double
SequenceTree::computeRobinsonFoulds(SequenceTree &t1, SequenceTree &t2){

  int numLeafs = t1.getNumLeafs();  
  if(numLeafs != t2.getNumLeafs()){
    USER_WARNING("trees have different num Leafs: " << numLeafs <<"!=" << t2.getNumLeafs());
    return -1;
  }

  SequenceTree::NodeVector nodes1;
  SequenceTree::NodeVector nodes2;
  
  //1. index map for leafs;
  str2int_hashmap name2index((int)(numLeafs*1.5));
  t1.addLeafs(nodes1);
  for(int i=0 ; i<numLeafs ; i++){
    name2index[NAME(nodes1[i])] = i;
  }

  //check that all leafs are the same in the two trees
 
  t2.addLeafs(nodes2);
  for(int i=0 ; i<numLeafs ; i++){
    str2int_hashmap::iterator iter = name2index.find(NAME(nodes2[i]));
    if(iter == name2index.end() ){
      USER_WARNING("trees have different leafs. " << NAME(nodes2[i]) <<" doesn't exist.");
      return -1;
    }
  }
  //2. compute splitt set for each tree.
  vector<BitVector> splitts1;
  t1.computeSplittSet(splitts1,nodes1,name2index);
  //t1.drawTree(cout);
  vector<BitVector> splitts2;
  t2.computeSplittSet(splitts2,nodes2,name2index);
  //t2.drawTree(cout);

  BitVectorPtr_set set1((int)(numLeafs*1.5));
  BitVectorPtr_set set2((int)(numLeafs*1.5));
    
 
  for(size_t i=0;i<splitts1.size();i++){
    if(!nodes1[i]->isLeaf() && !nodes1[i]->isRoot()){
      set1.insert(&splitts1[i]);
      //printSplitt(splitts1[i],name2index);
    }
  }
  for(size_t i=0;i<splitts2.size();i++){
    if(!nodes2[i]->isLeaf() && !nodes2[i]->isRoot()){
      set2.insert(&splitts2[i]);
      //printSplitt(splitts2[i],name2index);
    }
  }
  //3. compute RF formula over splitt sets;
  
  int in_1_notin_2=0;
  int in_2_notin_1=0;

  BitVectorPtr_set::iterator iter = set1.begin();
  for( ; iter!= set1.end() ; ++iter){
    if( set2.find(*iter) == set2.end()){
      //      PRINT("1 didn't find");PRINT(*iter);printSplitt(**iter,name2index);
      in_1_notin_2++;
    }
  }
  iter = set2.begin();
  for( ; iter!= set2.end() ; ++iter){
    if( set1.find(*iter) == set1.end()){
      //      PRINT("2 didn't find");PRINT(*iter); printSplitt(**iter,name2index);
      in_2_notin_1++;
    }
  }

  //PRINT(in_1_notin_2);PRINT(in_2_notin_1);

  //The Robinson Foulds distance between two trees is
  //the number of edges in one tree that are not in the other tree.
  //i.e RF(T1,T2) = |Splitts(T1)\Splitts(T2)| + |Splitts(T2)\Splitts(T1)|
  //The normalized meassure is thus RF(T1,T2)/(Splitts(T1)+Splitts(T2))
  //The number of splitts in a binary tree is (remember leaf edges are not splitts)
  //n-3, since the number of nodes in a rooted binary tree is n-1 and the root is not a
  //splitt node and also the two edges from the root describe the same splitt.

  return ((double) in_1_notin_2 + in_2_notin_1)/( set1.size() + set2.size());
 
}
float
SequenceTree::computeFloatRobinsonFoulds(SequenceTree &t1, SequenceTree &t2){

  int numLeafs = t1.getNumLeafs();
  if(numLeafs != t2.getNumLeafs()){
    USER_WARNING("trees have different num Leafs: " << numLeafs <<"!=" << t2.getNumLeafs());
    return -1;
  }

  SequenceTree::NodeVector nodes1;
  SequenceTree::NodeVector nodes2;

  //1. index map for leafs;
  str2int_hashmap name2index((int)(numLeafs*1.5));
  t1.addLeafs(nodes1);
  for(int i=0 ; i<numLeafs ; i++){
    name2index[NAME(nodes1[i])] = i;
  }

  //check that all leafs are the same in the two trees

  t2.addLeafs(nodes2);
  for(int i=0 ; i<numLeafs ; i++){
    str2int_hashmap::iterator iter = name2index.find(NAME(nodes2[i]));
    if(iter == name2index.end() ){
      USER_WARNING("trees have different leafs. " << NAME(nodes2[i]) <<" doesn't exist.");
      return -1;
    }
  }
  //2. compute splitt set for each tree.
  vector<BitVector> splitts1;
  t1.computeSplittSet(splitts1,nodes1,name2index);
  //t1.drawTree(cout);
  vector<BitVector> splitts2;
  t2.computeSplittSet(splitts2,nodes2,name2index);
  //t2.drawTree(cout);

  BitVectorPtr_set set1((int)(numLeafs*1.5));
  BitVectorPtr_set set2((int)(numLeafs*1.5));


  for(size_t i=0;i<splitts1.size();i++){
    if(!nodes1[i]->isLeaf() && !nodes1[i]->isRoot()){
      set1.insert(&splitts1[i]);
      //printSplitt(splitts1[i],name2index);
    }
  }
  for(size_t i=0;i<splitts2.size();i++){
    if(!nodes2[i]->isLeaf() && !nodes2[i]->isRoot()){
      set2.insert(&splitts2[i]);
      //printSplitt(splitts2[i],name2index);
    }
  }
  //3. compute RF formula over splitt sets;

  int in_1_notin_2=0;
  int in_2_notin_1=0;

  BitVectorPtr_set::iterator iter = set1.begin();
  for( ; iter!= set1.end() ; ++iter){
    if( set2.find(*iter) == set2.end()){
      //      PRINT("1 didn't find");PRINT(*iter);printSplitt(**iter,name2index);
      in_1_notin_2++;
    }
  }
  iter = set2.begin();
  for( ; iter!= set2.end() ; ++iter){
    if( set1.find(*iter) == set1.end()){
      //      PRINT("2 didn't find");PRINT(*iter); printSplitt(**iter,name2index);
      in_2_notin_1++;
    }
  }

  //PRINT(in_1_notin_2);PRINT(in_2_notin_1);

  //The Robinson Foulds distance between two trees is
  //the number of edges in one tree that are not in the other tree.
  //i.e RF(T1,T2) = |Splitts(T1)\Splitts(T2)| + |Splitts(T2)\Splitts(T1)|
  //The normalized meassure is thus RF(T1,T2)/(Splitts(T1)+Splitts(T2))
  //The number of splitts in a binary tree is (remember leaf edges are not splitts)
  //n-3, since the number of nodes in a rooted binary tree is n-1 and the root is not a
  //splitt node and also the two edges from the root describe the same splitt.

  return ((float) in_1_notin_2 + in_2_notin_1)/( set1.size() + set2.size());

}


void
SequenceTree::computeSplittSet(std::vector<BitVector> &splitts, SequenceTree::NodeVector &nodes,//ForEachSplitt,
		 str2int_hashmap &name2index){

  nodes.clear();
  recalcNodeIdsPostfixOrderAndAddInOrder(nodes);
  const size_t numNodes = nodes.size();
  const size_t numLeafs = getNumLeafs();
  splitts.clear();
  splitts.resize(numNodes);
  
  for(size_t i=0 ; i<numNodes ; i++){  
    if(nodes[i]->isLeaf()){
      int index = name2index[NAME(nodes[i])];
      splitts[i].setNumBits(numLeafs);
      splitts[i].setBit(index);
      splitts[i].flippAllIfPositionIsCleared(0);
    }
    else{
      SequenceTree::Node *child= nodes[i]->getRightMostChild();
      splitts[i] = splitts[child->getNodeId()];//set to first child
      child = child->getLeftSibling();
      for( ; child!=NULL; child = child->getLeftSibling() ){
	splitts[i].bitwiseEqual(splitts[child->getNodeId()]);//join the splitts
      }
    }
  }
  
}

void
SequenceTree::printSplitt(BitVector &splitt, str2int_hashmap &name2index, std::ostream& os){
  
  vector<string> index2name(name2index.size());
  
  str2int_hashmap::iterator iter = name2index.begin();
  for( ; iter!=name2index.end() ; ++iter){
    index2name[(*iter).second] = (*iter).first;
  }
  os << splitt << "  Hashcode: " << splitt.hashCode() << " "; 
  for(size_t i=0 ; i<splitt.getNumBits() ; i++){
    os << "(" << index2name[i] << ","<<splitt.getBit(i) << ") ";
  }
  os << endl;
}


//------------------------------------------------------------------------


void
SequenceTree::tree2distanceMatrix(StrDblMatrix &dm){
  int numLeafs = getNumLeafs();

  dm.resize(numLeafs);
  SequenceTree::NodeVector leafs;
  addLeafs(leafs);
  for(size_t i=0;i<leafs.size();i++){
    dm.setIdentifier(i,NAME(leafs[i]));
    for(size_t j=i+1;j<leafs.size();j++)
      dm.setDistance(i,j,0);
  }
  

  SequenceTree::NodeVector nodesOnPath;
  for(size_t i=0;i<leafs.size();i++){
    for(size_t j=i+1;j<leafs.size();j++){
      nodesOnPath.clear();
      addNodesOnPathExceptLCA(nodesOnPath,leafs[i],leafs[j]);
      double sum =0;
      for(size_t e=0;e<nodesOnPath.size();e++){
	sum+=EDGE(nodesOnPath[e]);
      }
      dm.setDistance(i,j,sum);
    }
  }
}
void
SequenceTree::tree2FloatdistanceMatrix(StrFloMatrix &fdm){
  int numLeafs = getNumLeafs();

  fdm.resize(numLeafs);
  SequenceTree::NodeVector leafs;
  addLeafs(leafs);
  for(size_t i=0;i<leafs.size();i++){
    fdm.setIdentifier(i,NAME(leafs[i]));
    for(size_t j=i+1;j<leafs.size();j++)
      fdm.setDistance(i,j,0);
  }


  SequenceTree::NodeVector nodesOnPath;
  for(size_t i=0;i<leafs.size();i++){
    for(size_t j=i+1;j<leafs.size();j++){
      nodesOnPath.clear();
      addNodesOnPathExceptLCA(nodesOnPath,leafs[i],leafs[j]);
      float sum =0;
      for(size_t e=0;e<nodesOnPath.size();e++){
	sum+=EDGE(nodesOnPath[e]);
      }
      fdm.setDistance(i,j,sum);
    }
  }
}


//--------------------------------------------------------------------

void
SequenceTree::createLeafNameToLeafIdMap(str2int_hashmap &name2id) const{
  SequenceTree::const_NodeVector leafs;
  addLeafs(leafs);
  name2id.clear();
  name2id.reserve((size_t)(1.5*leafs.size()));
  int leafId = 0;
  
  for(size_t i=0;i<leafs.size();i++)
    name2id[NAME(leafs[i])]=leafId++;
	    
}
//the leafs are added in order of their ids in name2id
void
SequenceTree::makeCanonical(const str2int_hashmap &name2id){
  SequenceTree::NodeVector leafs(getNumLeafs());
  
  SequenceTree::NodeVector  tmpvec;
  addLeafs(tmpvec);
  
  for(size_t i=0;i<tmpvec.size();i++){
    str2int_hashmap::const_iterator find = name2id.find(NAME(tmpvec[i]));
    if(find==name2id.end()){USER_WARNING("name doesn't exist: \"" << NAME(tmpvec[i])<<"\"");}
    leafs[(*find).second] = tmpvec[i];  
  }
  makeCanonical(leafs);
}










































