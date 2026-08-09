// Regression test for SequenceTree::mapSequencesOntoTree(std::istream&)
// (SequenceTree.cpp) - has no live caller anywhere in the codebase
// (its only two callers were dead simulator-related code, deleted -
// see project history), so this is the only way to verify this code
// path at all.
//
// One thing discovered writing this test, worth recording since
// nothing else in the project could have surfaced it without a live
// caller to exercise this code against: readSequenceLine() normalizes
// nucleotide characters to lowercase (via nucleotide2char()) - every
// expected string below is lowercase, not a copy-paste of the
// uppercase input.
//
// This reader used to be genuinely incompatible with the plain,
// unpadded continuation lines the rest of this project's PHYLIP files
// use (the ones Sequences2DistanceMatrix.cpp/Sequence.cpp already read
// correctly) - a peek-based name-field-skip blindly discarded 10 raw
// bytes whenever a continuation line didn't start with whitespace, on
// the assumption every line repeats a fixed-width name field (a
// stricter, older PHYLIP convention nothing in this project actually
// produces). Since this function has no live caller and that stricter
// convention added real incompatibility cost for no benefit, the skip
// was removed rather than kept - see SequenceTree.cpp's
// mapSequencesOntoTree(istream&). The multi-line fixtures below use
// plain, flush-left continuation lines (this project's actual
// convention) - the exact form that used to throw `Bad character`
// mid-read.

#undef NDEBUG
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include "fastphylo/core/SequenceTree.hpp"

using namespace std;

namespace
{

SequenceTree makeThreeLeafTree()
{
    istringstream newick("(Alpha:1,Beta:1,Gamma:1);");
    return SequenceTree(newick);
}

} // namespace

// Single-line sequences (no interleaving needed) - the baseline case.
static void test_single_line_sequences()
{
    SequenceTree tree = makeThreeLeafTree();

    istringstream phylip(" 3 8\n"
                         "Alpha     ACGTACGT\n"
                         "Beta      TTTTAAAA\n"
                         "Gamma     GGGGCCCC\n");
    tree.mapSequencesOntoTree(phylip);

    SequenceTree::NodeVector nodes;
    tree.addNodesInPrefixOrder(nodes);

    string alpha, beta, gamma;
    for (auto *n : nodes)
    {
        if (nodeName(n) == "Alpha")
        {
            alpha = nodeSeq(n);
        }
        if (nodeName(n) == "Beta")
        {
            beta = nodeSeq(n);
        }
        if (nodeName(n) == "Gamma")
        {
            gamma = nodeSeq(n);
        }
    }
    assert(alpha == "acgtacgt");
    assert(beta == "ttttaaaa");
    assert(gamma == "ggggcccc");
}

// Forces the interleaved-continuation loop (phylipReadInterleavedContinuation)
// to actually run, exercising the migrated shared implementation rather
// than just the opening-line read. Continuation lines are plain and
// flush-left, no leading whitespace at all - this project's actual
// convention, and the exact form that used to throw `Bad character`
// before the peek-based 10-char name-field skip was removed (see the
// file header comment).
static void test_multiline_interleaved_sequences()
{
    SequenceTree tree = makeThreeLeafTree();

    istringstream phylip(" 3 16\n"
                         "Alpha     ACGTACGTAC\n"
                         "Beta      TTTTAAAATT\n"
                         "Gamma     GGGGCCCCGG\n"
                         "ACGTAC\n"
                         "TTAATT\n"
                         "GGCCCC\n");
    tree.mapSequencesOntoTree(phylip);

    SequenceTree::NodeVector nodes;
    tree.addNodesInPrefixOrder(nodes);

    string alpha, beta, gamma;
    for (auto *n : nodes)
    {
        if (nodeName(n) == "Alpha")
        {
            alpha = nodeSeq(n);
        }
        if (nodeName(n) == "Beta")
        {
            beta = nodeSeq(n);
        }
        if (nodeName(n) == "Gamma")
        {
            gamma = nodeSeq(n);
        }
    }
    assert(alpha == "acgtacgtacacgtac");
    assert(beta == "ttttaaaattttaatt");
    assert(gamma == "ggggccccggggcccc");
}

// Two names in the file ("Fake1"/"Fake2") don't match any tree leaf,
// so both are redirected to the shared `garbage` buffer (Phase A of
// the consolidation plan). "Fake2"'s opening line is deliberately
// exactly `seqlen` characters (16) long - garbage.clear() before every
// name's read means garbage never accumulates across multiple
// unmatched names (each starts fresh), but a *single* unmatched line
// that happens to be exactly `seqlen` long would, if `garbage` weren't
// excluded from the completion check, make
// phylipAnySequenceComplete() falsely report done before any real
// sequence (still only 10/16 chars after the opening line) has
// actually finished - collapsing the continuation loop before it runs
// at all, leaving Alpha/Beta/Gamma short and triggering
// mapSequencesOntoTree(istream&)'s own trailing length-check error.
// This test fails loudly (an aborted process, via USER_ERROR) if that
// exclusion is ever lost.
static void test_unmatched_name_does_not_corrupt_real_sequences()
{
    SequenceTree tree = makeThreeLeafTree();

    istringstream phylip(" 5 16\n"
                         "Alpha     ACGTACGTAC\n"
                         "Beta      TTTTAAAATT\n"
                         "Fake1     AAAAAAAAAAAAAAAA\n"
                         "Gamma     GGGGCCCCGG\n"
                         "Fake2     CCCCCCCCCCCCCCCC\n"
                         "ACGTAC\n"
                         "TTAATT\n"
                         "AAAAAA\n"
                         "GGCCCC\n"
                         "CCCCCC\n");
    tree.mapSequencesOntoTree(phylip);

    SequenceTree::NodeVector nodes;
    tree.addNodesInPrefixOrder(nodes);

    string alpha, beta, gamma;
    for (auto *n : nodes)
    {
        if (nodeName(n) == "Alpha")
        {
            alpha = nodeSeq(n);
        }
        if (nodeName(n) == "Beta")
        {
            beta = nodeSeq(n);
        }
        if (nodeName(n) == "Gamma")
        {
            gamma = nodeSeq(n);
        }
    }
    assert(alpha == "acgtacgtacacgtac");
    assert(beta == "ttttaaaattttaatt");
    assert(gamma == "ggggccccggggcccc");
}

int main()
{ // NOLINT(bugprone-exception-escape) - a test binary crashing on an unexpected exception is a fine, loud failure mode.
    test_single_line_sequences();
    test_multiline_interleaved_sequences();
    test_unmatched_name_does_not_corrupt_real_sequences();
    cout << "SequenceTree_PhylipReader_test: all tests passed" << endl;
    return 0;
}
