//--------------------------------------------------
//
// File: SequenceTree.cpp
//
// Author: Mehmood Alam Khan,Isaac Elias
// e-mail: malagori@kth.se, isaac@nada.kth.se
//
//--------------------------------------------------

#include "fastphylo/core/SequenceTree.hpp"
#include "fastphylo/core/Exception.hpp"
#include "fastphylo/core/phylip_interleaved_reader.hpp"

#include <array>
#include <vector>
#include <string>

#include "fastphylo/core/InitAndPrintOn_utils.hpp"
#include <iostream>
#include "fastphylo/dna/string_compare.hpp"
#include "fastphylo/core/stl_utils.hpp"
#include "fastphylo/core/nucleotide.hpp"

using namespace std;

//=---------------------
void SequenceTree::printSequences(std::ostream &os)
{
    SequenceTree::NodeVector vec;
    vec.reserve(getNumNodes());
    addNodesInInfixOrder(vec);
    for (auto n : vec)
    {
        if (nodeName(n).empty())
        {
            continue;
        }
        os << n->data.s << endl;
    }
}

//---------------------
using str2node_map = std::unordered_map<const std::string, SequenceTree::Node *, hashstr, eqstr>;

void SequenceTree::mapSequencesOntoTree(std::vector<Sequence> &seqs)
{
    size_t numnodes = getNumNodes();
    SequenceTree::NodeVector nodes;
    nodes.reserve(numnodes);
    addNodesInPrefixOrder(nodes);

    // build hash map {name,node}
    str2node_map str2node(static_cast<int>(numnodes * 1.5));
    for (size_t i = 0; i < numnodes; i++)
    {
        if (!nodeName(nodes[i]).empty())
        {
            str2node[nodeName(nodes[i])] = nodes[i];
        }
    }

    // go through the nameseq pairs and look up node
    for (auto &seq : seqs)
    {
        auto iter = str2node.find(seq.name);
        if (iter != str2node.end())
        {
            nodeSeq((*iter).second).clear();
            nodeSeq((*iter).second).append(seq.seq);
        }
        else
        {
            //      USER_WARNING("Unknown name in tree: " << nameseqPairs[i]);
        }
    }
}

//--------------------
namespace
{
// Reads one line's worth of sequence characters from fin into *currseq,
// translating each via char2nucleotide() and stopping at end-of-line.
// Shared between mapSequencesOntoTree(istream&)'s first pass (each
// sequence's opening line) and its interleaved-continuation pass (later
// lines of the same, possibly multi-line, alignment) - both read one
// line into whichever sequence pointer is currently active. This was
// the main source of that function's cognitive complexity (61):
// duplicated inline, its nesting counted twice.
void readSequenceLine(std::istream &fin, string *currseq)
{
    while (true)
    {
        char c = fin.get();
        nucleotide n = char2nucleotide(c);
        if (DNA_NOT_ALLOWED == n)
        {
            if (isspace(c) == 0)
            {
                THROW_EXCEPTION("Bad character '" << c << "'");
            }
            else if (c != '\n')
            {
                continue; // skip space
            }
            break; // if '\n'
        }
        currseq->append(1, nucleotide2char(n));
    }
}

} // namespace

void SequenceTree::mapSequencesOntoTree(std::istream &fin)
{

    int numSequences;
    unsigned int seqlen;
    std::array<char, 100> tmp{};
    fin.getline(tmp.data(), 100);
    // NOLINTNEXTLINE(bugprone-unchecked-string-to-number-conversion) - checked below.
    if (sscanf(tmp.data(), "%d %d", &numSequences, &seqlen) != 2)
        THROW_EXCEPTION("expected \"<num sequences> <sequence length>\", got \"" << tmp.data() << "\"");

    // build hash map {name,node} and reserve mem for the sequences
    str2node_map str2node(static_cast<int>(getNumNodes() * 1.5));
    SequenceTree::NodeVector nodes;
    nodes.reserve(getNumNodes());
    addNodesInPrefixOrder(nodes);
    for (size_t i = 0; i < nodes.size(); i++)
    {
        nodeSeq(nodes[i]).clear();
        nodeSeq(nodes[i]).reserve(seqlen + 10);
        if (!nodeName(nodes[i]).empty())
        {
            str2node[nodeName(nodes[i])] = nodes[i];
        }
    }

    // names in the file are matched agains a node sequences If the name
    // doesn't exist in the tree then these sequences are read into a
    // garbage string.
    string garbage;
    garbage.reserve(seqlen + 10);
    std::vector<string> names(numSequences);
    std::vector<string *> sequences(numSequences, &garbage);

    // read the names and map the sequences onto the tree

    for (int i = 0; i < numSequences; i++)
    {
        fin >> names[i];
        garbage.clear();

        // find the node
        auto iter = str2node.find(string(names[i]));
        if (iter != str2node.end())
        {
            sequences[i] = &nodeSeq(((*iter).second));
        }

        readSequenceLine(fin, sequences[i]);
    }

    // read remaining sequences

    // The sequences aren't neccesarily on one line but my be spread out interleaving
    // over several lines. Therefore we read until seqlen chars have been read.
    phylipReadInterleavedContinuation(
        fin, numSequences, static_cast<size_t>(seqlen),
        [&sequences, &garbage](int i) -> size_t {
            // Only a real (matched) position may signal completion - one
            // still pointing at the shared `garbage` buffer must never
            // falsely trigger it, since every unmatched name in the file
            // appends into that same buffer, inflating its length faster
            // than any single real sequence grows (Phase A of
            // phylip_reader_consolidation_plan.md). Returning a value that
            // can never equal seqlen excludes it from the check.
            if (sequences[i] == &garbage)
            {
                return static_cast<size_t>(-1);
            }
            return sequences[i]->length();
        },
        [&](istream &in, int i) {
            // Skip any fully blank separator line(s) some interleaved PHYLIP
            // variants put between blocks - readSequenceLine() itself already
            // skips leading whitespace *within* a line generically, so no
            // further special-casing of the line's first character is
            // needed. This used to also blindly discard exactly 10 raw bytes
            // whenever a continuation line didn't start with whitespace, on
            // the assumption every continuation line repeats the (possibly
            // blank) fixed-width name field - a stricter, older PHYLIP
            // convention that made this reader unable to parse this
            // project's own plain, unpadded continuation lines (the ones
            // Sequences2DistanceMatrix.cpp/Sequence.cpp already read
            // correctly, and the only convention any of this project's real
            // fixtures actually use).
            char c = in.peek();
            while (c == '\n')
            {
                c = in.get();
            }

            garbage.clear();
            readSequenceLine(in, sequences[i]);
        });

    // CHECK THAT ALL STRINGS HAVE THE SAME LENGTH
    for (size_t i = 0; i < nodes.size(); i++)
    {
        if (nodeSeq(nodes[i]).length() != seqlen && !nodeSeq(nodes[i]).empty())
        {
            THROW_EXCEPTION("Sequence not of correct length: " << nodes[i]->data.s.name);
        }
    }
}

//--------------------------------------------------------------------

// the leafs are added in order of their ids in name2id
void SequenceTree::makeCanonical(const str2int_hashmap &name2id)
{
    SequenceTree::NodeVector leafs(getNumLeafs());
    SequenceTree::NodeVector tmpvec;
    addLeafs(tmpvec);

    for (size_t i = 0; i < tmpvec.size(); i++)
    {
        auto find = name2id.find(nodeName(tmpvec[i]));
        if (find == name2id.end())
        {
            std::cerr << "warning: name doesn't exist: \"" << nodeName(tmpvec[i]) << "\"" << std::endl;
        }
        leafs[(*find).second] = tmpvec[i];
    }
    makeCanonical(leafs);
}
