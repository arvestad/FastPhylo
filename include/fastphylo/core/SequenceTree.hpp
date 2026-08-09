//--------------------------------------------------
//
// File: SequenceTree.hpp
//
// Author: Mehmood Alam Khan, Isaac Elias
// e-mail: malagori@kth.se, isaac@nada.kth.se
//
//--------------------------------------------------
#pragma once

#include "fastphylo/core/Tree.hpp"
#include <string>
#include "fastphylo/core/InitAndPrintOn_utils.hpp"
#include "fastphylo/core/Object.hpp"
#include "fastphylo/core/Sequence.hpp"

//---------------------------------------
// A REGULAR SEQUENCE TREE
// Each node has a struct of a sequence and a double.
//--------------------------------------

// Accessors for a SequenceTree node's payload (a Sequence_double: a
// Sequence + an edge length). Free function templates, not member
// methods on TreeNode<Data,...> itself - Data varies across the
// generic Tree<> framework's other instantiations, so these only make
// sense for the specific Data shape SequenceTree uses. Kept as
// templates (not hardcoded to SequenceTree::Node*) since
// NeighborJoining.hpp's own templates call these while still generic
// over the node type, matching the old macros' purely structural
// (duck-typed) behavior exactly. Return by reference so existing
// assignment call sites (e.g. NAME(n) = "...") keep working unchanged.
template <class Node> inline auto &nodeName(Node *node)
{
    return node->data.s.name;
}
template <class Node> inline auto &nodeSeq(Node *node)
{
    return node->data.s.seq;
}
template <class Node> inline auto &nodeEdge(Node *node)
{
    return node->data.dbl;
}

class SequenceTree : public Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double>>
{
  public:
    // the same constructers as in Tree
    SequenceTree() : Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double>>()
    {
    }

    SequenceTree(Sequence_double d)
        : Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double>>(d)
    {
    }

    SequenceTree(const SequenceTree &t)
        : Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double>>(t)
    {
    }

    SequenceTree(const TreeNode<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double>> &n)
        : Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double>>(n)
    {
    }

    SequenceTree(char *newickstr)
        : Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double>>(newickstr)
    {
    }
    SequenceTree(std::istream &in)
        : Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double>>(in)
    {
    }

    template <class Data2, class DataInit2, class DataPrintOn2>
    SequenceTree(const Tree<Data2, DataInit2, DataPrintOn2> &t, Sequence_double defaultData)
        : Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double>>(t, defaultData)
    {
    }

    virtual ~SequenceTree()
    {
    }

    //----------------------------------
    // ADDITION PRINT FUNCTIONS
    // The print functions in Tree are not shadowed.
    void printSequences(std::ostream &os);

    //---------------------------
    // MAKE CANONICAL
    // the leafs are added in order of their ids in name2id, returns false
    // if the name2id map doesn't contain all leafs
    void makeCanonical(const str2int_hashmap &name2id);
    // the makeCanonical in the super class
    using Tree<Sequence_double, Data_init<Sequence_double>, Data_printOn<Sequence_double>>::makeCanonical;

    //--------
    void mapSequencesOntoTree(std::vector<Sequence> &seqs);

    // Takes a phylip sequence file and maps the sequences onto the tree
    // nodes by looking up the name.
    // It is possible to do the mapping from several files. That is some names
    // may exist in one file and other names in another file.
    void mapSequencesOntoTree(std::istream &in);
};

using tree2int_map = std::unordered_map<const SequenceTree, int, objhash, objeq>;
