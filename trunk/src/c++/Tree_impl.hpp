//--------------------------------------------------
//                                        
// File: Tree_impl.hpp                             
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: Tree_impl.hpp,v 1.50 2006/12/19 12:36:23 isaac Exp $                                 
//
// Implementation of Tree template in Tree.hpp 
//--------------------------------------------------
#ifndef TREE_IMPL_HPP
#define TREE_IMPL_HPP

#include <string>
#include <iostream>
#include <sstream>
#include <istream>
#include <stack>
#include <assert.h>
#include "log_utils.hpp"
#include <algorithm>
#include "xml_output_global.hpp"

#define TREE_ASSERT(stmnt) if(true){stmnt;}

//-------------------------------------------------------
// CONSTRUCTORS


TREE_TEMPLATE TREE::Tree(char *newickstr){
  _nullVariables();
  std::string str(newickstr);
  std::istringstream in(str);
  root = this->initSubtreeFromStream(in);
  TREE_ASSERT(assertTreeStructure());
}

TREE_TEMPLATE TREE::Tree(Data d){
  _nullVariables();
  root = new TREENODE(d,this);
  recalcNodeStructure();
  TREE_ASSERT(assertTreeStructure());
}

TREE_TEMPLATE TREE::Tree(const TREE &t){
   _nullVariables();
   (*this) = t;
  TREE_ASSERT(assertTreeStructure());
}

//creates a copy of the subtree
TREE_TEMPLATE TREE::Tree(const TREENODE &n){
  root = new TREENODE(n);
  numNodes = 100;
  recalcNodeStructure();
  TREE_ASSERT(assertTreeStructure());

}


TREE_TEMPLATE template<class Data2, class DataInit2, class DataPrintOn2>
TREE::Tree(const Tree<Data2,DataInit2,DataPrintOn2> &t, Data defaultData){
  root = new TREENODE(*(t.getRoot()),defaultData);
  numNodes = t.getNumNodes();
  recalcNodeStructure();
  TREE_ASSERT(assertTreeStructure());
}

TREE_TEMPLATE TREE& TREE::operator=(const TREE &t){
  if ( root != NULL )
    delete root;
  if( t.root == NULL )
    root = NULL;
  else{
  root = new TREENODE(*(t.root));
  
  nodeId = t.nodeId;

  NodeVector vec;
  vec.reserve(numNodes);
  addNodesInPrefixOrder(vec);
  const_NodeVector vec2;
  vec2.reserve(numNodes);
  t.addNodesInPrefixOrder(vec2);
  for ( size_t i = 0 ; i < vec.size() ; i++ ){
    TREENODE *n = vec[i];
    n->ownertree = this;
    n->nodeId = vec2[i]->nodeId;
  }

  }
  numNodes = t.numNodes;
  numLeafs = t.numLeafs;


  TREE_ASSERT(assertTreeStructure());

  return *this;
}

TREE_TEMPLATE TREE::~Tree(){
  if ( root != NULL )
    delete root;
}

  
//---------------------------------------------------------
// READING FROM STREAM

TREE_TEMPLATE std::istream&
TREE::objInitFromStream(std::istream &is){
  if ( root != NULL )
    delete root;
  _nullVariables();

  root = initSubtreeFromStream(is);
  
  return is;
}


TREE_TEMPLATE TREENODE *
TREE::initSubtreeFromStream(std::istream &in){
  //  cout << "parse: " << peek(f) << endl;
  
  //the newick string format of the node is either
  //1. "A" followed by "," or ")"  the name may be empty but then there can't be an edge length
  //2. "A:17" followed by "," or ")"
  //3. "(...):17 followed by "," or ")"
  //GENERAL FORMAT
  //Every node is either just DATA i.e. a leaf or (CHILDREN)DATA i.e. an inner node.
  //The data after the inner node may not be separated by whitespace from the ')'. If it is
  //this has to be handled by the DataInitializer.
  
  skipWhiteSpace(in);
  
  //LEAF 
   if ( in.peek() != '(' ) {
     Data d;
     dataInit(in,d);

     numLeafs++;
     numNodes++;
     return new TREENODE(d,this);
        
   }
   //INNER NODE
   else{
     in.get();//skip the ´(´
     skipWhiteSpace(in);

     if ( in.peek() == ')' ){
       std::cerr << "WRONG INPUT FOR TREE" << std::endl;
       exit(1);
     }

     //pointer to the rightmost child
     TREENODE *rightChild = initSubtreeFromStream(in);
     skipWhiteSpace(in);

     while ( in.peek() != ')' ){
       if ( in.peek() != ',' ){
         std::cerr << "WRONG INPUT FOR TREE" << std::endl;
         exit(1);
       }
       in.get();//skip the ','
      
       TREENODE *tmp = initSubtreeFromStream(in);
       rightChild->rightSibling = tmp;
       tmp->leftSibling = rightChild;
       rightChild = tmp;
      
       skipWhiteSpace(in);
     }

    
     in.get();//skip the ')'
     skipWhiteSpace(in);
     Data d;
     dataInit(in,d);
     TREENODE *node = new TREENODE(d,this);
     numNodes++;
     node->rightMostChild = rightChild;
     while ( rightChild != NULL ){
       rightChild->parent = node;
       rightChild = rightChild->leftSibling;
     }
     return node;
   }
}


//---------------------------------------------------------
// PRINTING AND DRAWING
TREE_TEMPLATE std::ostream&
TREE::printOn(std::ostream &os) const{

  if(xmlPrint) {
    if(root != NULL)
         os << root;      
  } else
    {
      if(root==NULL)
	return os << "NULLTREE";

      return os << root <<";";
    }
 

}


TREE_TEMPLATE void
TREE::drawSubtree(std::ostream &os, const TREENODE *n,std::string &prefix) const{
  if ( n == root ){
    if ( root == NULL ){
      os << "NULL TREE\n";
      return;
    }
    os << "[" << n->nodeId << "] ";
    dataPrintOn(os,n->data) << "\n";
    prefix.append("*");
    const TREENODE *c = n->getRightMostChild();
    while ( c != NULL ){
      this->drawSubtree(os, c, prefix);
      if ( c->getLeftSibling() != NULL )
        os << prefix << "\n";
      c = c->getLeftSibling();
    }
  }
  else if ( n->isLeaf() ){
    os << prefix << "-------["<< n->nodeId<< "] ";
    dataPrintOn(os,n->data) << "\n";
  }
  else{
    int len = prefix.size();
  
    os << prefix << "-------["<< n->nodeId<< "] ";
    dataPrintOn(os, n->data) << "\n";
    if ( n->getLeftSibling() == NULL && prefix.size() > 0){
      prefix.replace(prefix.size()-1,sizeof(char),1,' ');
    }
    prefix.append("       *");
    const TREENODE *c = n->getRightMostChild();
    while ( c != NULL ){
      this->drawSubtree(os, c, prefix);
      if ( c->getLeftSibling() != NULL )
        os << prefix << "\n";
      c = c->getLeftSibling();
    }

    prefix.erase(len);
  }
}


TREE_TEMPLATE void
TREE::drawSubtreeIDs(std::ostream &os, const TREENODE *n,std::string &prefix) const{
  if ( n == root ){
    if ( root == NULL ){
      os << "NULL TREE\n";
      return;
    }
    os << "[" << n->nodeId << "]\n";
    prefix.append("*");
    const TREENODE *c = n->getRightMostChild();
    while ( c != NULL ){
      this->drawSubtreeIDs(os, c, prefix);
      if ( c->getLeftSibling() != NULL )
        os << prefix << "\n";
      c = c->getLeftSibling();
    }
  }
  else if ( n->isLeaf() ){
    os << prefix << "-------["<< n->nodeId<< "]\n";
  }
  else{
    int len = prefix.size();
  
    os << prefix << "-------["<< n->nodeId<< "]\n";
    if ( n->getLeftSibling() == NULL && prefix.size() > 0){
      prefix.replace(prefix.size()-1,sizeof(char),1,' ');
    }
    prefix.append("       *");
    const TREENODE *c = n->getRightMostChild();
    while ( c != NULL ){
      this->drawSubtreeIDs(os, c, prefix);
      if ( c->getLeftSibling() != NULL )
        os << prefix << "\n";
      c = c->getLeftSibling();
    }

    prefix.erase(len);
  }
}

TREE_TEMPLATE void TREE::recalcNodeStructure(){

  nodeId = 0;
  int numVisited = 0;
  int numLeafsVisited = 0;

  NodeVector vec;
  vec.reserve(numNodes);
  addNodesInPrefixOrder(vec);
  for ( size_t i = 0 ; i < vec.size() ; i++ ){
    TREENODE *n = vec[i];
    n->ownertree = this;
    n->nodeId = getNewNodeId();

    numVisited++;
    if ( n->isLeaf()) numLeafsVisited++; 
  }

  numNodes = numVisited;
  numLeafs = numLeafsVisited;
}


TREE_TEMPLATE bool
TREE::isBinary() const {

  const_NodeVector vec;
  vec.reserve(numNodes);
  addNodesInPrefixOrder(vec);
  for ( size_t i = 0 ; i < vec.size() ; i++ ){
    int deg = vec[i]->getDegree();
    if ( deg != 1 && deg != 3 )
      return false;
  }

  return true;
}


//--------------------------------------------------------------
// CHANGING TREE STRUCTURE

TREE_TEMPLATE void 
TREE::moveNode(TREENODE *node, 
	       TREENODE *newparent){
  if( node->isRoot() ){
    USER_ERROR("Node is root and cannot be moved");
    return;
  }
  if( newparent->isDescendantOf(node) ){
    USER_ERROR("Node cannot be moved to one of its children");
    return;
  } 
  node->detachFromParent();
  node->parent = newparent;
  node->rightSibling = NULL;

  if( newparent->isLeaf() ){
    newparent->rightMostChild = node;
    node->leftSibling = NULL;
  } else {
    newparent->rightMostChild->rightSibling = node;
    node->leftSibling = newparent->rightMostChild;
    newparent->rightMostChild = node;
  }  
  recalcNodeStructure();
  TREE_ASSERT(assertTreeStructure());
}


TREE_TEMPLATE void TREE::joinTreeAtRoot(TREE &t){
  joinTreeAtNode(root,t);
}

TREE_TEMPLATE void
TREE::joinTreeAtNode(TREENODE *n, TREE &t){

  
  TREENODE *newnode = t.root;
  newnode->parent = n;
  newnode->rightSibling = NULL;
  
  
  if ( n->isLeaf() ){
    n->rightMostChild = newnode;
    newnode->leftSibling = NULL;
  }
  else{
    n->rightMostChild->rightSibling = newnode;
    newnode->leftSibling = n->rightMostChild;
    n->rightMostChild = newnode;
  }
  
  recalcNodeStructure();

  
  t.numNodes = 0;
  t.numLeafs = 0;
  t.root = NULL;
  TREE_ASSERT(t.assertTreeStructure());
  TREE_ASSERT(assertTreeStructure());
}

TREE_TEMPLATE void
TREE::joinTreeAtNodes(TREENODE *n,
                      TREENODE *node_other_tree){
  node_other_tree->ownertree->reRootAt(node_other_tree);
  joinTreeAtNode(n,*(node_other_tree->ownertree));

  TREE_ASSERT(assertTreeStructure());
}

TREE_TEMPLATE TREENODE *
TREE::insertNodeOnPathToParent(TREENODE *n, Data newparentdata){

  if ( n->isRoot() )
    USER_ERROR("can insert node on path to parent for the root" );

  TREENODE *parent = n->parent;
  //remove n the parent.
  n->detachFromParent();

  //----
  //create the parent node with default data.
  TREENODE *newparent = new TREENODE(newparentdata,this);
  newparent->parent = parent;
  newparent->nodeId = getNewNodeId();
  numNodes++;
  newparent->rightSibling = NULL;
  if ( parent->rightMostChild != NULL )
    parent->rightMostChild->rightSibling = newparent;
  newparent->leftSibling = parent->rightMostChild;
  parent->rightMostChild = newparent;
  //----

  //add n as child to the new parent
  newparent->rightMostChild = n;
  n->leftSibling = NULL;
  n->parent = newparent;

  TREE_ASSERT(assertTreeStructure());
  return newparent;  
}

TREE_TEMPLATE void
TREE::splittTree(TREENODE *new_tree_root, TREE &new_tree){
  assert( new_tree_root->getTree() == this );
  if ( new_tree.root != NULL )
    delete new_tree.root;
  
  if ( new_tree_root == root ){
    root = NULL;
    numNodes = 0;
    numLeafs = 0;
  }
  else{
    new_tree_root->detachFromParent();
    recalcNodeStructure();
  }
  new_tree.root = new_tree_root;
  new_tree.recalcNodeStructure();
  TREE_ASSERT(assertTreeStructure());
  TREE_ASSERT(new_tree.assertTreeStructure());
}

TREE_TEMPLATE TREENODE*
TREE::detachFromParentAndAddAsSiblings(TREENODE *a, TREENODE *b, Data newparentdata){

  assert ( a->parent == b->parent );

  TREENODE *parent = a->parent;
  //remove a and b from the children of their parent.
  a->detachFromParent();
  b->detachFromParent();
  //----
  //create the parent node with default data.
  TREENODE *newparent = new TREENODE(newparentdata,this);
  newparent->parent = parent;
  newparent->nodeId = getNewNodeId();
  numNodes++;
  newparent->rightSibling = NULL;
  if ( parent->rightMostChild != NULL )
    parent->rightMostChild->rightSibling = newparent;
  newparent->leftSibling = parent->rightMostChild;
  parent->rightMostChild = newparent;
  //----

  //add a and b as children to the new parent
  newparent->rightMostChild = a;
  a->leftSibling = b;
  b->rightSibling = a;
  a->parent = newparent;
  b->parent = newparent;

  TREE_ASSERT(assertTreeStructure());
  return newparent;
}

TREE_TEMPLATE void
TREENODE::detachFromParent(){

  assert ( !isRoot() );
  
  if ( isRightMostChild() )
    parent->rightMostChild = leftSibling;
  else
    rightSibling->leftSibling = leftSibling;
    
  if ( leftSibling != NULL )
    leftSibling->rightSibling = rightSibling;
 
  leftSibling = NULL;
  rightSibling = NULL;
  parent = NULL;
}


TREE_TEMPLATE void
TREE::removeAndDelete(TREENODE *n){

  if ( n->isRoot() ){
    root = NULL;
    delete n;
    numNodes = 0;
    numLeafs = 0;
  }
  else{
    if ( n->isRightMostChild() )
      n->parent->rightMostChild = n->leftSibling;
    else
      n->rightSibling->leftSibling = n->leftSibling;
    
    if ( n->leftSibling != NULL )
      n->leftSibling->rightSibling = n->rightSibling;
        
    delete n;
    recalcNodeStructure();
  }

  TREE_ASSERT(assertTreeStructure());
}



TREE_TEMPLATE void
TREE::reRootAt(TREENODE *n){
  if ( n->ownertree != this ){
    PROG_ERROR("not correct owner tree");
  }
  if ( n == root )
    return;
  //WARNING the root is NOT considered to be a leaf 
  if ( n->isLeaf() )
    numLeafs--;
  if ( root->rightMostChild->leftSibling == NULL )
    numLeafs++;
  
  n->setAsRoot();
  root = n;
  TREE_ASSERT(assertTreeStructure());
}

TREE_TEMPLATE void
TREENODE::setAsRoot(){
  if ( isRoot() )
    return;

  //1. Ask parent to reroot
  parent->setAsRoot();

  //2. Remove from parents children
  if ( isRightMostChild() )
    parent->rightMostChild = leftSibling;
  else
    rightSibling->leftSibling = leftSibling;
    
  if ( leftSibling != NULL )
    leftSibling->rightSibling = rightSibling;

  //3. Add parent to the children
  if ( ! isLeaf() )
    rightMostChild->rightSibling = parent;

  parent->leftSibling = rightMostChild;

  rightMostChild = parent;
  parent->parent = this;
  
  //4. Fix so that the node doesn't have any siblings or parent
  leftSibling = NULL;
  rightSibling = NULL;
  parent = NULL;
}

TREE_TEMPLATE bool
TREE::reRootAtHighDegreeNode(){

  NodeVector vec;
  addNodesInPrefixOrder(vec);
  size_t i = 0;
  for  (  ; i < vec.size() ; i++ ){
    if ( vec[i]->getDegree() > 1 ){
      reRootAt(vec[i]);
      return true;
    }
  }

  //if there are only two nodes in the tree
  return false;
}

TREE_TEMPLATE TREENODE *
TREE::shortcutNode(TREENODE *n){
  if ( n->isLeaf() ){
    TREENODE *p = n->parent;
    removeAndDelete(n);
    return p->rightMostChild;
  }
  if ( n->isRoot() ){
    root = n->rightMostChild;
    TREENODE *oldrmchild = root->rightMostChild;
    root->parent = NULL;
    //if the child is not an only child
    if ( root->leftSibling != NULL ){
      TREENODE *leftMostSibling = n->getLeftMostChild();
      if ( oldrmchild != NULL )
      oldrmchild->rightSibling = leftMostSibling;
      else numLeafs--;
      leftMostSibling->leftSibling = oldrmchild;
      root->rightMostChild = root->leftSibling;
      root->leftSibling->rightSibling = NULL;
      root->leftSibling = NULL;
    }
    //set all the parent pointers
    TREENODE *newchild = root->rightMostChild;
    for ( ; newchild != oldrmchild ; newchild = newchild->leftSibling )
      newchild->parent = root;
    n->rightMostChild = NULL;
    delete n;
    numNodes--;
    TREE_ASSERT(assertTreeStructure());
    return oldrmchild;
  }
  else{//not the root
    TREENODE *parent = n->parent;
    n->detachFromParent();
    TREENODE *oldrmchild = parent->rightMostChild;
    //if n is not the only child of the parent
    if ( oldrmchild != NULL ){
      TREENODE *leftMostSibling = n->getLeftMostChild();
      oldrmchild->rightSibling = leftMostSibling;
      leftMostSibling->leftSibling = oldrmchild;
    }
    parent->rightMostChild = n->rightMostChild;
    //set all the parent pointers
    TREENODE *newchild = parent->rightMostChild;
    for ( ; newchild != oldrmchild ; newchild = newchild->leftSibling )
      newchild->parent = parent;
    n->parent = NULL;
    n->rightMostChild = NULL;
    n->leftSibling = NULL;
    n->rightSibling = NULL;
    delete n;
    numNodes--;  
    TREE_ASSERT(assertTreeStructure());
    return oldrmchild;
  }

}

 //------
// COLLAPSE
//remove the node and add its children to the parent.
TREE_TEMPLATE void
TREE::collapse(TREENODE *n){

  //if leaf nothing is done
    if ( getNumNodes() == 1 || n->getDegree()==1 ){
      return;
    }
    if ( n->parent == NULL ){
      PRINT(n->parent==NULL);
      reRootAt(n->getRightMostChild());
    }

    TREENODE *parent = n->parent;
    n->detachFromParent();
    TREENODE *left= n->rightMostChild;
    for ( ; left->leftSibling != NULL ; left = left->leftSibling )
      left->parent = parent;
    left->parent = parent;
    left->leftSibling = parent->rightMostChild;
    parent->rightMostChild->rightSibling = left;
    parent->rightMostChild = n->rightMostChild;
    numNodes--;
    n->rightMostChild = NULL;
    delete n;
    
    TREE_ASSERT(assertTreeStructure());
}

TREE_TEMPLATE void
TREE::collapseChildren(TREENODE *n){
    NodeVector vec;
    TREENODE *c = n->rightMostChild;
    for ( ; c != NULL ; c = c->leftSibling )
      vec.push_back(c);
    for ( size_t i=0 ; i<vec.size() ; i++ )
      collapse(vec[i]);
    TREE_ASSERT(assertTreeStructure());
    
}


//------------------------------------------------------------------
TREE_TEMPLATE bool TREE::assertTreeStructure() const{

 
  const_NodeVector prefixvec;
  addNodesInPrefixOrder(prefixvec);

  const_NodeVector postfixvec;
  addNodesInPrefixOrder(postfixvec);

  ASSERT_EQ( postfixvec.size() , prefixvec.size());
  ASSERT_EQ( postfixvec.size() , numNodes );

  size_t numLeafsVisited = 0;
 
  size_t i = 0 ;

  for ( ; i < numNodes ; i++ ){
    const TREENODE *n = prefixvec[i];
    ASSERT_EQ ( (void *) n->ownertree, (void*)this);
    if ( n->rightSibling != NULL )
      ASSERT_EQ ( (void *) n->parent,(void *) n->rightSibling->parent );
    if ( !n->isLeaf() )
      assert ( n->rightMostChild->rightSibling == NULL );

    if ( n->isLeaf()) numLeafsVisited++;
  }
  
  ASSERT_EQ ( numLeafsVisited , numLeafs );

  return true;
}


//--------------------------------------------
// MAKE CANONICAL
TREE_TEMPLATE
struct nodeidcmp
{
  bool operator()(const TREENODE *const &n1, const TREENODE * const& n2) const
  {
    return (n1->getNodeId()>n2->getNodeId());
  }
};

TREE_TEMPLATE void 
TREE::makeCanonical(const std::vector<TREENODE *> &leafs){
    
  if(leafs.size()==0) return;

  //clear the nodeId counter and assign ids to leafs
  nodeId = 0;
  for(size_t i=0;i<leafs.size();i++){
    if(((void*)leafs[i]->getTree())!=((void*)this))
      PROG_ERROR("leafs not from this tree");

    leafs[i]->nodeId = getNewNodeId();
  }
  
  //reroot at parent of leaf 0
  TREENODE *newroot = NULL;
  if( leafs[0]->getParent()==NULL ){
    TREENODE *newroot = leafs[0]->getLeftMostChild();
    if( newroot==NULL )
      return;//there is only one node in the whole tree
  }
  else
    newroot = leafs[0]->getParent();

  reRootAt(newroot);

  //
  NodeVector nodes;
  nodes.reserve(getNumNodes());
  addNodesInPostfixOrder(nodes);
  NodeVector children; 
  nodeidcmp<Data_type,DataInitializer_type, DataPrintOn_type> cmp;
    
  for(size_t i=0 ; i<nodes.size() ; i++){
    if(nodes[i]->isLeaf()) continue;

    nodes[i]->nodeId = getNewNodeId();

    //sort children according to their node ids.
    children.clear();  
    nodes[i]->addChildren(children);
    std::sort(children.begin(),children.end(),cmp);

    nodes[i]->rightMostChild = children[0];
    TREENODE *prevSibling = children[0];
    children[0]->rightSibling = NULL;
    for(size_t c=1 ; c<children.size() ; c++){
      children[c]->rightSibling = prevSibling;
      prevSibling->leftSibling = children[c];
      prevSibling = children[c];
    }
    children[children.size()-1]->leftSibling = NULL;
  }


  TREE_ASSERT(assertTreeStructure());
}

// EQUALS
TREE_TEMPLATE bool 
TREE::equals(const Object *o) const{
  if ( getNumNodes() != ((const TREE *)o)->getNumNodes() )
    return false;
  
  const_NodeVector nodes1;
  nodes1.reserve(getNumNodes());
  addNodesInPrefixOrder(nodes1);
  const_NodeVector nodes2;
  nodes2.reserve(getNumNodes());
  ((const TREE *)o)->addNodesInPrefixOrder(nodes2);
  
  for(size_t i=0 ; i<nodes1.size() ; i++){
    if( nodes1[i]->getNodeId()!=nodes2[i]->getNodeId() )
      return false;
  }

  return true;
}

// HASHCODE
TREE_TEMPLATE size_t 
TREE::hashCode() const{
  const_NodeVector nodes;
  nodes.reserve(getNumNodes());
  addNodesInPrefixOrder(nodes);

  size_t hash = 5381;
  
  for(size_t i=0 ; i<nodes.size() ; i++ )
    hash = ((hash << 5) + hash) + nodes[i]->getNodeId(); /* hash * 33 + c */

  //  PRINT(hash);drawTreeIDs(std::cout);
  
  return hash;
}

//-------------------------------------------------------
// TreeNode Implementation
//-------------------------------------------------------


TREE_TEMPLATE TREENODE::TreeNode(Data d, TREE *owner) : data(d){
  ownertree = owner;
  nodeId = owner->getNewNodeId();
  parent = NULL;
  rightSibling = NULL;
  leftSibling = NULL;
  rightMostChild = NULL;
}

TREE_TEMPLATE TREENODE::TreeNode(const TREENODE &n) : data(n.data) {

  if ( n.isLeaf() ){
    parent = NULL;
    rightMostChild = NULL;
    rightSibling = NULL;
    leftSibling = NULL;
  }
  else{
    parent = NULL;
    leftSibling = NULL;
    rightSibling = NULL;
    
    TREENODE *child = n.rightMostChild;
    TREENODE *cpy = new TREENODE(*child);

    cpy->parent = this;
    cpy->rightSibling = NULL;
    rightMostChild = cpy;

    child = child->leftSibling;
    while ( child != NULL ){
      TREENODE *nextcpy = new TREENODE(*child);
      nextcpy->parent = this;
      cpy->leftSibling = nextcpy;
      nextcpy->rightSibling = cpy;
      cpy = nextcpy;
      child = child->leftSibling;
    }
    cpy->leftSibling = NULL;
  }
}



TREE_TEMPLATE TREENODE::~TreeNode(){
  if ( ! isLeaf() ){
    TREENODE *child = rightMostChild;
    while ( child != NULL ){
      TREENODE *tmp = child;
      child = child->leftSibling;
      delete tmp;
    }
  }  
}

TREE_TEMPLATE template<class Data2, class DataInit2, class DataPrintOn2>
TREENODE::TreeNode(const TreeNode<Data2,DataInit2,DataPrintOn2> &n, Data defaultData) : data(defaultData) {

  if ( n.isLeaf() ){
    rightMostChild = NULL;
    rightSibling = NULL;
    leftSibling = NULL;
  }
  else{
    parent = NULL;
    leftSibling = NULL;
    rightSibling = NULL;
    
    const TreeNode<Data2,DataInit2,DataPrintOn2> *child = n.getRightMostChild();
    TREENODE *cpy = new TREENODE(*child,defaultData);

    cpy->parent = this;
    cpy->rightSibling = NULL;
    rightMostChild = cpy;

    child = child->getLeftSibling();
    while ( child != NULL ){
      TREENODE *nextcpy = new TREENODE(*child,defaultData);
      nextcpy->parent = this;
      cpy->leftSibling = nextcpy;
      nextcpy->rightSibling = cpy;
      cpy = nextcpy;
      child = child->getLeftSibling();
    }
    cpy->leftSibling = NULL;
  }
  
}


TREE_TEMPLATE int
TREENODE::getNumChildren() const{

  int numC = 0;
  const TREENODE *c = getRightMostChild();
  for ( ; c!=NULL ; c = c->leftSibling )
    numC++;

  return numC;
}
//-------------------
TREE_TEMPLATE size_t
TREENODE::getDegree() const{
  size_t degree = 0;
  if ( parent != NULL )
    degree++;
  TREENODE *child = rightMostChild;
  for ( ; child != NULL ; child = child->leftSibling )
    degree++;

  return degree;
}

TREE_TEMPLATE std::ostream &
TREENODE::printOn(std::ostream &os) const{


  if ( isLeaf() ){
    if(xmlPrint) {
      os << "<leaf";
    }
    std::ostringstream outstr;
    ownertree->dataPrintOn(outstr, data);
    std::string str = outstr.str();
    if ( str.size() == 0 )
      str = std::string("n") + nodeId;
    if(xmlPrint) {
      return os << str << "</leaf>";
    } else
    {
    return os << str;
    }

  }
  else{
    //PENDING SLOW

    const TREENODE *child = rightMostChild;



    if(xmlPrint) {
      os << "<branch" ;
       ownertree->dataPrintOn(os, data); 
       os << child;



    } 
    else   {
      os <<"(";
      os << child;
    }



    child = child->leftSibling;
    
    for ( ; child != NULL ; child = child->leftSibling )
      {
	if(xmlPrint) {
	  os <<  child;
	} else
	  {
	    os << "," << child;
	  }

      }

	if(xmlPrint) {
          os << "</branch>";
	} 
	else {
          os << ")";
	}
        

	  if(!xmlPrint) {  ownertree->dataPrintOn(os, data); }

        
    return os;
  }
  
}

TREE_TEMPLATE TREENODE *
TREENODE::addChild(Data d){
  TREENODE *node = new TREENODE(d,ownertree);
  node->parent = this;
  node->rightMostChild = NULL;
  node->rightSibling = NULL;

  if ( ! this->isLeaf() )
    ownertree->numLeafs++;
  ownertree->numNodes++;

  if ( rightMostChild != NULL )
  rightMostChild->rightSibling = node;
  node->leftSibling = rightMostChild;
  rightMostChild = node;

  TREE_ASSERT(ownertree->assertTreeStructure());
  return node;
}


//------------------------------------
TREE_TEMPLATE TREENODE *
TREENODE::getRightMostSibling(){
  TREENODE *tmp = rightSibling;

  if ( tmp == NULL )
    return NULL;

  while ( tmp->rightSibling != NULL )
    tmp = tmp->rightSibling;

  return tmp;
}

TREE_TEMPLATE const TREENODE *
TREENODE::getRightMostSibling() const{
  const TREENODE *tmp = rightSibling;

  if ( tmp == NULL )
    return NULL;

  while ( tmp->rightSibling != NULL )
    tmp = tmp->rightSibling;

  return tmp;
}
//------------------------------------
TREE_TEMPLATE const TREENODE *
TREENODE::getLeftMostSibling() const{
  const TREENODE *tmp = leftSibling;

  if ( tmp == NULL )
    return NULL;

  while ( tmp->leftSibling != NULL )
    tmp = tmp->leftSibling;

  return tmp;
}

TREE_TEMPLATE TREENODE *
TREENODE::getLeftMostSibling(){
  TREENODE *tmp = leftSibling;

  if ( tmp == NULL )
    return NULL;

  while ( tmp->leftSibling != NULL )
    tmp = tmp->leftSibling;

  return tmp;
}
//------------------------------------



TREE_TEMPLATE const TREENODE *
TREENODE::getLeftMostChild() const {
  const TREENODE *tmp = rightMostChild;
  if ( tmp == NULL )
    return NULL;

  while ( tmp->leftSibling != NULL )
    tmp = tmp->leftSibling;
  

  return tmp;
}


TREE_TEMPLATE TREENODE *
TREENODE::getLeftMostChild() {
  TREENODE *tmp = rightMostChild;
  if ( tmp == NULL )
    return NULL;

  while ( tmp->leftSibling != NULL )
    tmp = tmp->leftSibling;
  

  return tmp;
}

//-------------------------------------------
TREE_TEMPLATE const TREENODE *
TREENODE::getRightSiblingOfAncestor() const {

  const TREENODE *tmp = this;
  while ( tmp->rightSibling == NULL ){
    tmp = tmp->parent;
    if ( tmp == NULL )
      return NULL;
  }

  return tmp->rightSibling;
}


TREE_TEMPLATE TREENODE *
TREENODE::getRightSiblingOfAncestor(){

  TREENODE *tmp = this;
  while ( tmp->rightSibling == NULL ){
    tmp = tmp->parent;
    if ( tmp == NULL )
      return NULL;
  }

  return tmp->rightSibling;
}

//----------------------------------------

TREE_TEMPLATE const TREENODE *
TREENODE::getRightMostDescendantLeaf() const{

  const TREENODE *tmp = this;
  while ( tmp->rightMostChild != NULL )
    tmp = tmp->rightMostChild;

  return tmp;
}

TREE_TEMPLATE TREENODE *
TREENODE::getRightMostDescendantLeaf(){

  TREENODE *tmp = this;
  while ( tmp->rightMostChild != NULL )
    tmp = tmp->rightMostChild;

  return tmp;
}

//----------------------------------------
TREE_TEMPLATE const TREENODE *
TREENODE::getLeftMostDescendantLeaf() const{

  const TREENODE *tmp = this;
  while ( ! tmp->isLeaf() )
    tmp = tmp->getLeftMostChild;

  return tmp;
}

TREE_TEMPLATE TREENODE *
TREENODE::getLeftMostDescendantLeaf(){

  TREENODE *tmp = this;
  while ( ! tmp->isLeaf() )
    tmp = tmp->getLeftMostChild;

  return tmp;
}

//------------------------------------------------------------------------------
//-------------------------------------------------
// RETRIVING THE NODES IN PARITCULAR ORDER.
// 
TREE_TEMPLATE void
TREE::addNodesInPrefixOrder(std::vector<TREENODE *> &nodes, TREENODE *n) {
  if ( n == NULL ) return;
  
  nodes.push_back(n);

  if ( !n->isLeaf() ){
    TREENODE *child = n->rightMostChild;
    for ( ; child != NULL ; child = child->leftSibling )
      addNodesInPrefixOrder(nodes,child);    
  }   
}


TREE_TEMPLATE void
TREE::addNodesInPrefixOrder(std::vector<const TREENODE *> &nodes, const TREENODE *n) {
  if ( n == NULL ) return;
  
  nodes.push_back(n);

  if ( !n->isLeaf() ){
    const TREENODE *child = n->rightMostChild;
    for ( ; child != NULL ; child = child->leftSibling )
      addNodesInPrefixOrder(nodes,child);    
  }   
}

TREE_TEMPLATE void
TREE::addNodesInPostfixOrder(std::vector<TREENODE *> &nodes, TREENODE *n){
  if ( n == NULL ) return;
  
  if ( !n->isLeaf() ){
    TREENODE *child = n->rightMostChild;
    for ( ; child != NULL ; child = child->leftSibling )
      addNodesInPostfixOrder(nodes,child);    
  }
  
  nodes.push_back(n);
}

TREE_TEMPLATE void
TREE::addNodesInPostfixOrder(std::vector<const TREENODE *> &nodes, const TREENODE *n) {
  if ( n == NULL ) return;
  
  if ( !n->isLeaf() ){
    const TREENODE *child = n->rightMostChild;
    for ( ; child != NULL ; child = child->leftSibling )
      addNodesInPostfixOrder(nodes,child);    
  }
  
  nodes.push_back(n);
}

TREE_TEMPLATE void
TREE::addNodesInPrefixOrderLeftRight(std::vector<TREENODE *> &nodes, TREENODE *n) {
  if ( n == NULL ) return;
  
  nodes.push_back(n);

  if ( !n->isLeaf() ){
    TREENODE *child = n->getLeftMostChild();
    for ( ; child != NULL ; child = child->rightSibling )
      addNodesInPrefixOrderLeftRight(nodes,child);    
  }   
}

TREE_TEMPLATE void
TREE::addNodesInInfixOrder(std::vector<TREENODE *> &nodes, TREENODE *n){
  if ( n == NULL ) return;
  
  if ( n->isLeaf() ){
    nodes.push_back(n);
  }
  else {
    TREENODE *child = n->rightMostChild;
    addNodesInInfixOrder(nodes,child);

    nodes.push_back(n);
    
    child = child->leftSibling;
    for ( ; child != NULL ; child = child->leftSibling )
      addNodesInInfixOrder(nodes,child);    
  }
}

TREE_TEMPLATE void
TREE::addNodesInInfixOrder(std::vector<const TREENODE *> &nodes, const TREENODE *n){
  if ( n == NULL ) return;
  
  if ( n->isLeaf() ){
    nodes.push_back(n);
  }
  else {
    const TREENODE *child = n->rightMostChild;
    addNodesInInfixOrder(nodes,child);

    nodes.push_back(n);
    
    child = child->leftSibling;
    for ( ; child != NULL ; child = child->leftSibling )
      addNodesInInfixOrder(nodes,child);    
  }
}


TREE_TEMPLATE void
TREE::addLeafs(std::vector<TREENODE *> &nodes, TREENODE *n){
 if ( n == NULL ) return;
   
  if ( n->isLeaf() )
    nodes.push_back(n);
  else{
    TREENODE *child = n->rightMostChild;
    for ( ; child != NULL ; child = child->leftSibling )
      addLeafs(nodes,child);    
  }
}

TREE_TEMPLATE void
TREE::addLeafs(std::vector<const TREENODE *> &nodes, const TREENODE *n){
 if ( n == NULL ) return;
   
  if ( n->isLeaf() )
    nodes.push_back(n);
  else{
    const TREENODE *child = n->rightMostChild;
    for ( ; child != NULL ; child = child->leftSibling )
      addLeafs(nodes,child);    
  }
}


TREE_TEMPLATE void
TREE::addInternalNodes(std::vector<TREENODE *> &nodes, TREENODE *n){
  if ( n == NULL ) return;
  
  if ( !n->isLeaf() ){
    nodes.push_back(n);
    TREENODE *child = n->rightMostChild;
    for ( ; child != NULL ; child = child->leftSibling )
      addInternalNodes(nodes,child);    
  }
}


TREE_TEMPLATE void
TREE::addNodesWithDegree(std::vector<TREENODE *> &nodes, TREENODE *n, size_t degree){
  if ( n == NULL ) return;
  
  if ( n->getDegree() == degree )
    nodes.push_back(n);
  
  if ( ! n->isLeaf() ){
    TREENODE *child = n->rightMostChild;
    for ( ; child != NULL ; child = child->leftSibling )
      addNodesWithDegree(nodes,child,degree);    
  }
}


TREE_TEMPLATE void
TREE::addNodesOnPath(std::vector<TREENODE *> &nodes, TREENODE *n1, TREENODE *n2){

  NodeVector tmpvec1;
  TREENODE *tmp1 = n1;
  do {
    tmpvec1.push_back(tmp1);
    tmp1 = tmp1->parent;
  } while ( tmp1 != NULL );
  
  NodeVector tmpvec2;
  TREENODE *tmp2 = n2;
  do {
    tmpvec2.push_back(tmp2);
    tmp2 = tmp2->parent;
  } while ( tmp2 != NULL );

  int i = tmpvec1.size()-1;
  int j = tmpvec2.size()-1;

  while ( tmpvec1[i] == tmpvec2[j] ){
      i--;
      j--;
  }

  for ( int k = 0 ; k <= i+1 ; k++ )
    nodes.push_back(tmpvec1[k]);

  while ( j >= 0 )
    nodes.push_back(tmpvec2[j--]);
  
}



TREE_TEMPLATE void
TREE::addNodesOnPathExceptLCA(std::vector<TREENODE *> &nodes, TREENODE *n1, TREENODE *n2){

  NodeVector tmpvec1;
  TREENODE *tmp1 = n1;
  do {
    tmpvec1.push_back(tmp1);
    tmp1 = tmp1->parent;
  } while ( tmp1 != NULL );
  
  NodeVector tmpvec2;
  TREENODE *tmp2 = n2;
  do {
    tmpvec2.push_back(tmp2);
    tmp2 = tmp2->parent;
  } while ( tmp2 != NULL );

  int i = tmpvec1.size()-1;
  int j = tmpvec2.size()-1;

  while ( tmpvec1[i] == tmpvec2[j] ){
      i--;
      j--;
  }

  for ( int k = 0 ; k <= i ; k++ )
    nodes.push_back(tmpvec1[k]);

  while ( j >= 0 )
    nodes.push_back(tmpvec2[j--]);
  
}






TREE_TEMPLATE void
TREE::addNodesInIdOrder(std::vector<const TREENODE *> &nodes) const{
  nodes.resize(nodeId+1,NULL);//some positions in the vector may remain NULL if the id is unused
  for(size_t i=0;i<nodes.size();i++)
    nodes[i]=NULL;
  
  const_NodeVector tmpnodes;
  addNodesInPrefixOrder(tmpnodes);
  for(size_t i=0;i<tmpnodes.size();i++)
    nodes[tmpnodes[i]->nodeId] =tmpnodes[i];
}









// //
// // THE ALGORITHM 
// //1. assign leaf ids [0,numleafs]
// //2. for i=0...numnodes
// //      if node with id i only has one neighbor assign it the next id

// TREE_TEMPLATE void
// TREE::assignNodeIdsInUniqueWayFromLeafs(const std::vector<TREENODE *> &leafs){
  
//   //clear all node ids
//   NodeVector allnodes;
//   addNodesInPostfixOrder(allnodes);
//   for(size_t i=0;i<allnodes.size();i++)
//     allnodes[i]->nodeId = -1;


//   //the id counter starts at 1 so the first element is empty  
//   NodeVector nodesInIdOrder(getNumNodes()+1);
//   for(size_t i=0;i<nodesInIdOrder.size();i++)
//     nodesInIdOrder[i]=NULL;

//   //clear the nodeId counter and assign ids to leafs
//   nodeId = 0;
//   for(size_t i=0;i<leafs.size();i++){
//     if(((void*)leafs[i]->getTree())!=((void*)this))
//       PROG_ERROR("leafs not from this tree");

//     leafs[i]->nodeId = getNewNodeId();
//     nodesInIdOrder[i] = leafs[i];
//     //PRINT(i);PRINT(leafs[i]->nodeId);
//   }
  
  
//   //do the assignment of the internal nodes according to the algorithm
//   int currentId = 0;
//   NodeVector neighs;
//   while(currentId<getNumNodes()){
//     assert(currentId<nodeId);//make sure that the id we are looking for isn't bigger than the biggest assigned
//     //get the unassigned neighbor
//     assert(nodesInIdOrder[currentId]!=NULL);
//     neighs.clear();
//     nodesInIdOrder[currentId]->addNeighbors(neighs);

//     TREENODE *unassignedNeighbor =NULL;
//     int numUnassignedNeighs = 0;
//     for(size_t i=0;i<neighs.size();i++)
//       if(neighs[i]->nodeId == -1 ){
// 	unassignedNeighbor = neighs[i];
// 	numUnassignedNeighs++;
//       }
//     //if only has one neighbor give it the next id
//     if(numUnassignedNeighs == 1){
//       unassignedNeighbor->nodeId = getNewNodeId();
//       nodesInIdOrder[unassignedNeighbor->nodeId ] = unassignedNeighbor;
//     }

//     currentId++;//go to next node
//   }
  
//   TREE_ASSERT(assertTreeStructure());
  
// }



// TREE_TEMPLATE void
// TREE::getEulerianWalk(std::vector<int> &idsonpath) const{
  
//   const_NodeVector nodes;
//   addNodesInIdOrder(nodes);
//   const_NodeVector tmp;
//   //sortedNeighbors[i] contains a vector of all neighbors of node i
//   //sorted in decreasing node id order.
//   std::vector<const_NodeVector> sortedNeighbors(nodes.size(),tmp);
//   nodeidcmp<Data_type,DataInitializer_type, DataPrintOn_type> cmp;

//   for(size_t i=0;i<nodes.size();i++){
//     if(nodes[i]==NULL)continue;
//     nodes[i]->addNeighbors(sortedNeighbors[i]);
//     std::sort(sortedNeighbors[i].begin(),sortedNeighbors[i].end(),cmp);
//   }

//   std::stack<const TREENODE *> nodeStack;
//   std::vector<bool> isInStack(nodes.size(),false);
//   const TREENODE *current = nodes[0];
//   //SEPARATOR();
//   while(true){
//     idsonpath.push_back(current->nodeId);
//     //PRINT(current->nodeId);
//     isInStack[current->nodeId] = true;
//     const TREENODE *nextCurrent = NULL;
//     while(!sortedNeighbors[current->nodeId].empty()){
//       const TREENODE *tmp = sortedNeighbors[current->nodeId].back();
//       sortedNeighbors[current->nodeId].pop_back();
//       if(!isInStack[tmp->nodeId]){
// 	nextCurrent = tmp;
// 	break;
//       }
//     }
//     if(nextCurrent==NULL){
//       if(nodeStack.empty()) break;
//       current = nodeStack.top();
//       nodeStack.pop();
//     } else{
//       nodeStack.push(current);
//       current = nextCurrent;
//     }
//   }
// }

// TREE_TEMPLATE size_t 
// TREE::getHashCodeFromEulearianWalk(const std::vector<int> &idsonpath){
//   size_t hash = 5381;

//   std::vector<int>::const_iterator iter = idsonpath.begin();
//   for(;iter!=idsonpath.end();++iter)
//     hash = ((hash << 5) + hash) + *iter; /* hash * 33 + c */
  
//   return hash;
// }

// TREE_TEMPLATE bool 
// TREE::eulerianWalksEquals(const std::vector<int> &p1, const std::vector<int> &p2){
//   if(p1.size()!=p2.size())return false;

//   for(size_t i=0;i<p1.size();i++)
//     if(p1[i]!=p2[i])
//       return false;

//   return true;
// }







#endif // TREE_IMPL_HPP

















