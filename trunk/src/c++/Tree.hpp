//--------------------------------------------------
//                                        
// File: Tree.hpp                              
//                             
// Author: Isaac Elias         
// e-mail: isaac@nada.kth.se   
//                             
// cvs: $Id: Tree.hpp,v 1.26 2006/12/19 09:24:36 isaac Exp $                                 
//
// A generic Tree where each node has userdefined data.
//--------------------------------------------------


#ifndef TREE_HH
#define TREE_HH

#include <string>
#include <iostream>
#include <sstream>
#include <istream>
#include <vector>

#include "Object.hpp"
#include "stl_utils.hpp"
#include "file_utils.hpp"
#include "log_utils.hpp"
#include "InitAndPrintOn_utils.hpp"

//
// The Tree class has a reference to a TreeNode representing the root
// of the tree. Each TreeNode contains a pointer to its parent, right
// most child, and the left and right sibling.
// 
// Each TreeNode contains one Data template with the user specified
// data. In addition to this parameter, to declare a tree a
// DataInitializer and DataPrintOn function has to be supplied. These
// functions are used for respectively reading and writing the Data to
// a stream.  Some typical print and read functions are declared in
// InitAndPrintOn_utils.hpp. The functions may write and read any
// character except for ',' ')', '(', and ';'.
//
// Example:
// Tree<int, Data_init<int>, Data_printOn<int> >
// Tree<XXX> // is the same as Tree<XXX, empty_Data_init<XXX>, empty_Data_printOn<XXX> >
//

#define TREE Tree<Data,DataInitializer,DataPrintOn>
#define TREENODE TreeNode<Data,DataInitializer,DataPrintOn>
#define TREE_TEMPLATE template<class Data, class DataInitializer, class DataPrintOn>


//full decleration at end of this file.
TREE_TEMPLATE class TreeNode;


template<class Data, class DataInitializer = empty_Data_init<Data>, class DataPrintOn = empty_Data_printOn<Data> >
class Tree : public Object
{
public:
  //types used in the tree
  typedef TREENODE Node;
  typedef std::vector<const TREENODE *> const_NodeVector;
  typedef std::vector<TREENODE *> NodeVector;
  typedef Data Data_type;
  typedef DataInitializer DataInitializer_type;
  typedef DataPrintOn DataPrintOn_type;

  //---------------
  //CONSTRUCTORS
  Tree(){_nullVariables();}
  //creates a tree with only the root
  Tree(Data d);
  Tree(const TREE &t);
  //creates a copy of the subtree
  Tree(const TREENODE &n);

  
  //GENERAL FORMAT
  //Ex. "(leaf1,leaf2)parent" 
  //Every node is either just DATA i.e. a leaf or (CHILDREN)DATA i.e. an inner node.
  //The CHILDREN nodes are seperated by ','.
  //The data after the inner node may not be separated by whitespace from the ')'. If it is
  //this has to be handled by the DataInitializer.
  Tree(char *newickstr);
  Tree(std::istream &in){_nullVariables(); root = initSubtreeFromStream(in);}

  virtual ~Tree();
  TREE& operator=(const TREE&t);
  

  //Init the tree from a stream.
  //Inherited from Object. 
  virtual std::istream& objInitFromStream(std::istream &is);
  
  //-----
  //Copy structure of another type of template tree.
  //The data is not initialized.
  template<class Data2, class DataInit2, class DataPrintOn2>
  Tree(const Tree<Data2,DataInit2,DataPrintOn2> &t, Data defaultData);
  

  //----
  //THE ROOT
  //Even if the root only has one neighbor it is not considered to be a leaf.
  //Therefore if you want the tree to work as a regular tree call reRootAtHighDegreeNode()
  TREENODE *getRoot(){ return root;}
  const TREENODE *getRoot()const { return root;}
  //reroot the tree
  void reRootAt(TREENODE *n);
  //returns false if there are only two nodes
  bool reRootAtHighDegreeNode();

  
  //The nodes a and b are detached from their parent and then added
  //to ancestor as siblings with a parent node inserted. (As in NJ)
  //Their new parent node is returned.
  TREENODE* detachFromParentAndAddAsSiblings(TREENODE *a, TREENODE *b, Data newparentdata);
                                
  //removes the subtree rooted at n and deletes it. The data in the
  //nodes is still owned by the caller.
  void removeAndDelete(TREENODE *n);
  
 
  int getNumLeafs() const { return numLeafs;}
  int getNumNodes() const { return numNodes;}
  
  //returns true if tree is binary
  bool isBinary() const;
  
  //--------------------------------
  // PRINTING
  
  // prints in the formate
  std::ostream& printOn(std::ostream& os) const;

  //------------------
  //DRAWING
  // draws tree structure with each node id 
  void drawTreeIDs(std::ostream& os)const{ std::string prefix; this->drawSubtreeIDs(os,root,prefix);}
  void drawSubtreeIDs(std::ostream& os, const TREENODE *node, std::string &prefix) const;
  
  void drawTree(std::ostream& os) const{ this->drawSubtree(os,root);}
  void drawSubtree(std::ostream& os, const TREENODE *node, std::string &prefix) const;
  void drawSubtree(std::ostream& os, const TREENODE *node) const{
    std::string prefix; this->drawSubtree(os,node,prefix);
  };

  //-------------------------------------------------------
  // CHANGING THE NODE STRUCTURE

  //-------
  // MERGING TREES
  //The nodes of the other tree is transfered into the current tree.
  void joinTreeAtRoot(TREE &t);
  void joinTreeAtNode(TREENODE *n, TREE &t);
  void joinTreeAtNodes(TREENODE *node_this_tree, TREENODE *node_other_tree);


  //------- 
  //Detaches node from its parent and adds it as a child to
  //newparent, which is another node in the tree.  Doesn't work if
  //node is the root.
  void moveNode(TREENODE *node, TREENODE *newparent);

  //-----
  //SPLITT OFF A SUBTREE
  //The input tree new_tree is assigned the new_tree_root as root. 
  void splittTree(TREENODE *new_tree_root, TREE &new_tree);
  
  //Inserts a node between n and the parent of n.
  TREENODE *insertNodeOnPathToParent(TREENODE *n, Data newparentdata);
  //----------------
  //SHORT CUT
  //n is removed from the tree and deleted. Its children are added
  //as children to the parent of n. The rightmost child of n is added as rightmost child
  // to the parent of n and the old children of the parent comes after the children of n.
  //If n is the root the rightmost child of n becomes the new root and the other children are
  //added to this node.
  //The old rightmost child of the parent is returned (this simplifies iteration over the new children).
  TREENODE * shortcutNode(TREENODE *n);

  //------
  // COLLAPSE
  //remove the node and add its children to the parent.
  //If the node is a leaf nothing happens
  void collapse(TREENODE *n);
  void collapseChildren(TREENODE *n);
  

  //-------------------------------------------------------
  // RETRIEVE NODES IN SPECIFIC ORDERS.
  // Adds the nodes in the subtree to the vector

  //PREFIX ORDER
  static void addNodesInPrefixOrder(std::vector<TREENODE *> &nodes, TREENODE *n);
  void addNodesInPrefixOrder(std::vector<TREENODE *> &nodes){ addNodesInPrefixOrder(nodes,root);}
  static void addNodesInPrefixOrder(std::vector<const TREENODE *> &nodes, const TREENODE *n);
  void addNodesInPrefixOrder(std::vector<const TREENODE *> &nodes)const{ addNodesInPrefixOrder(nodes,root);}
  //A prefix traversals in which the left children are traversed before the right
  static void addNodesInPrefixOrderLeftRight(std::vector<TREENODE *> &nodes, TREENODE *n);
  void addNodesInPrefixOrderLeftRight(std::vector<TREENODE *> &nodes){ addNodesInPrefixOrderLeftRight(nodes,root);}
  
  //POSTFIX ORDER
  static void addNodesInPostfixOrder(std::vector<TREENODE *> &nodes, TREENODE *n);
  void addNodesInPostfixOrder(std::vector<TREENODE *> &nodes){ addNodesInPostfixOrder(nodes,root);}
  static void addNodesInPostfixOrder(std::vector<const TREENODE *> &nodes, const TREENODE *n) ;
  void addNodesInPostfixOrder(std::vector<const TREENODE *> &nodes)const{ addNodesInPostfixOrder(nodes,root);}

  //INFIX ORDER
  static  void addNodesInInfixOrder(std::vector<TREENODE *> &nodes, TREENODE *n);
  void addNodesInInfixOrder(std::vector<TREENODE *> &nodes){ addNodesInInfixOrder(nodes,root);}
  static  void addNodesInInfixOrder(std::vector<const TREENODE *> &nodes, const TREENODE *n);
  void addNodesInInfixOrder(std::vector<const TREENODE *> &nodes)const { addNodesInInfixOrder(nodes,root);}

  // LEAFS
  //the root is not a leaf even if it only has one neighbor
  static void addLeafs(std::vector<TREENODE *> &nodes, TREENODE *n);
  void addLeafs(std::vector<TREENODE *> &nodes){ addLeafs(nodes,root);};
  static void addLeafs(std::vector<const TREENODE *> &nodes, const TREENODE *n);
  void addLeafs(std::vector<const TREENODE *> &nodes) const { addLeafs(nodes,root);};
  
  //INTERNAL NODES (pending const functions)
  static void addInternalNodes(std::vector<TREENODE *> &nodes, TREENODE *n);
  void addInternalNodes(std::vector<TREENODE *> &nodes){ addInternalNodes(nodes,root);};

  //NODES OF SPECIFIC DEGREE (pending const functions)
  static void addNodesWithDegree(std::vector<TREENODE *> &nodes, TREENODE *n, size_t degree);
  void addNodesWithDegree(std::vector<TREENODE *> &nodes, size_t degree){ addNodesWithDegree(nodes,root,degree);};

  //NODES ON PATH (pending const functions)
  void addNodesOnPath(std::vector<TREENODE *> &nodes, TREENODE *n1, TREENODE *n2);
  void addNodesOnPathExceptLCA(std::vector<TREENODE *> &nodes, TREENODE *n1, TREENODE *n2);


  //------------------------------------------------
  //RECOMPUTE THE IDS OF THE TREE
  void recalcNodeStructure();
  void recalcNodeIdsPostfixOrderAndAddInOrder(std::vector<TREENODE *> &nodes){
    nodes.clear();
    addNodesInPostfixOrder(nodes);
    nodeId = 0;
    for(size_t i=0 ; i<nodes.size() ; i++)
      nodes[i]->nodeId = nodeId++;
  }


  //some positions in the vector may remain NULL if the id is unused
  void addNodesInIdOrder(std::vector<const TREENODE *> &nodes) const;

  //-------- 

  
  //---------------------------------------------------------
  // MAKE CANONICAL
  // The leafs are given ids according to their order in the node vector.
  // The tree is then rerooted at the parent of leaf with id 0.
  // Subsequently, internal nodes are given the id of the descendant leaf with lowest id.
  //
  void makeCanonical(const std::vector<TREENODE *> &leafs);
  
  // EQUALS
  // Checks if the datastructures are identical with respect to the root
  // and node ids. To check if two trees have the same topology the first call
  // makeCanonical on both trees.
  virtual bool equals(const Object *o) const;
  // HASH CODE 
  virtual size_t hashCode() const;


  //just for checking that there are no bugs
  bool assertTreeStructure() const;

  //--------------------- PRIVATE ---------------------------
private:
  friend class TREENODE;
  //Used for printing and writing
  DataInitializer dataInit;
  DataPrintOn dataPrintOn;

  //the root node. Also unrooted trees can be used. The root node is just
  //the place from which the traversal iterators start.
  TREENODE *root;

  //a counter for giving each node in the tree a unique id.
  int nodeId;
  int getNewNodeId(){ return nodeId++;};
  
  size_t numNodes;
  size_t numLeafs;

  //Nulls all the instance variables.
  void _nullVariables(){
    root = NULL; nodeId = 0;
    numNodes = 0; numLeafs = 0;
  }

  
  //reads a subtree from the stream and the root node is returned.
  TREENODE *initSubtreeFromStream(std::istream &in);


};



//-------------------------------------------------------
// TreeNode declaration
//-------------------------------------------------------



TREE_TEMPLATE
class TreeNode : public Object{
public:  
  typedef Data Data_type;
  typedef DataInitializer DataInitializer_type;
  typedef DataPrintOn DataPrintOn_type;
  typedef TREE Tree_type;

  
  //the user data
  Data data;

  int getNodeId() const {return nodeId;}
  TREE *getTree() const {return ownertree;}

  //----------- SPECIFIC NODE QUESTIONS ------------------------  
  //the root is not a leaf nomatter if it only has one neighbor
  bool isLeaf() const { return rightMostChild == NULL;}
  bool isRoot() const { return parent == NULL;}

  bool isRightMostChild() const { return ( !isRoot() && parent->rightMostChild == this);}
  bool isDescendantOf(TREENODE *n){
    if(getParent()==NULL) return false;
    else if(getParent()==n) return true;
    else return getParent()->isDescendantOf(n);
  }

  int getNumChildren() const;
  //degree
  size_t getDegree() const;
  
  // Create a child of this node whose data is d.
  TREENODE *addChild(Data d);
  
  //printing a subtree.
  std::ostream& printOn(std::ostream& os) const;

  //reRoot
  void setAsRoot();
  
    
  //----------- TRAVERSALS OF NODE POINTERS ------------------------
  const TREENODE *getParent() const {return parent;}
  TREENODE *getParent() {return parent;}

  const TREENODE *getLeftSibling() const {return leftSibling;}
  TREENODE *getLeftSibling() {return leftSibling;};
  
  const TREENODE *getRightSibling() const {return rightSibling;}  
  TREENODE *getRightSibling() {return rightSibling;};

  const TREENODE *getLeftMostSibling() const;
  TREENODE *getLeftMostSibling();
  
  const TREENODE *getRightMostSibling() const;
  TREENODE *getRightMostSibling();
  
  const TREENODE *getRightMostChild() const {return rightMostChild;}
  TREENODE *getRightMostChild() {return rightMostChild;};
  
  const TREENODE *getLeftMostChild() const;
  TREENODE *getLeftMostChild();

  //finds the node immediately to the right on the same level or above.
  const TREENODE *getRightSiblingOfAncestor() const;
  TREENODE *getRightSiblingOfAncestor();
  
  const TREENODE *getRightMostDescendantLeaf() const;
  TREENODE *getRightMostDescendantLeaf();
  
  const TREENODE *getLeftMostDescendantLeaf() const;
  TREENODE *getLeftMostDescendantLeaf();

  //GETTING NEIGHBOURS AND THE CHILDREN INTO NODEVECTOR
  void addNeighbors(std::vector<const TREENODE *> &nodes) const {
    if(parent!=NULL) nodes.push_back(parent);
    addChildren(nodes);
  }
  void addNeighbors(std::vector<TREENODE *> &nodes) {
    if(parent!=NULL) nodes.push_back(parent);
    addChildren(nodes);
  }

  void addChildren(std::vector<const TREENODE *> &nodes) const{
    
    const TREENODE *c=rightMostChild;
    while(c!=NULL) {
      nodes.push_back(c);
      c = c->leftSibling;
    }
  }
  
  void addChildren(std::vector<TREENODE *> &nodes) {
    
    TREENODE *c=rightMostChild;
    while(c!=NULL) {
      nodes.push_back(c);
      c = c->leftSibling;
    }
  }
  

 
  //-------------- PRIVATE ---------------------
private:  
  
  friend class TREE;

  TREE *ownertree;
 
  TREENODE *parent;
  TREENODE *leftSibling;
  TREENODE *rightSibling;
  TREENODE *rightMostChild;
  
  //every node in a tree is given a unique id
  int nodeId;

  //---
  // CONSTRUCTORS
  // These are private since there needs to be a owner tree.
  TreeNode(Data d, TREE *ownertree);
  //creates a copy of the subtree the owner tree is not set
  TreeNode(const TREENODE &n);
  TREENODE& operator=(const TREENODE &n){ PROG_ERROR("Not implemented"); return *this;}

  //deletes all children
  virtual ~TreeNode();

  //copy node structure of other template type
  //creates a copy of the subtree the owner tree is not set
  template<class Data2, class DataInit2, class DataPrintOn2>
  TreeNode(const TreeNode<Data2,DataInit2,DataPrintOn2> &n, Data defaultData);

  
  //removes the node from the parent list. The node can not be the root.
  void detachFromParent();

};


//---------------------
// INCLUDE THE IMPLEMENTATION
#include "Tree_impl.hpp"
//----------------------------


#endif // TREE_HH







