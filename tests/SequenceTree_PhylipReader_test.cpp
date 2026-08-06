// Regression test for SequenceTree::mapSequencesOntoTree(std::istream&)
// (SequenceTree.cpp), the one of the three legacy interleaved-PHYLIP
// readers with no live caller in any shipped binary - see
// phylip_reader_consolidation_plan.md. Written as a real, committed,
// CMake-wired test as part of that plan's Phase C, since this is the
// only way to verify this code path at all (its only other caller,
// Simulator.cpp, is itself only reachable from the unbuilt
// simulated_phylogenies/ tools).
//
// Two things discovered writing this test, worth recording since
// nothing else in the project could have surfaced them without a live
// caller to exercise this code against:
// - readSequenceLine() normalizes nucleotide characters to lowercase
//   (via nucleotide2char()) - every expected string below is
//   lowercase, not a copy-paste of the uppercase input.
// - This reader's continuation-line convention is genuinely
//   incompatible with the plain, unpadded continuation lines the rest
//   of this project's PHYLIP files use (the ones
//   Sequences2DistanceMatrix.cpp/Sequence.cpp read correctly) - its
//   peek-based name-field-skip logic (Phase A's difference #4) treats
//   any continuation line starting with non-whitespace as having a
//   stray 10-character name field to discard, so a plain
//   already-verified-working example file, fed to this reader
//   instead, throws mid-read. Confirmed by direct experiment - a
//   plain continuation line throws `Bad character`; the same content
//   padded with 10 leading spaces (as if every continuation line
//   repeated a blank name field) reads correctly. All multi-line
//   fixtures below use that space-padded form, since it's the only
//   one this reader actually supports.

#undef NDEBUG
#include <cassert>
#include <iostream>
#include <sstream>
#include <string>
#include "fastphylo/core/SequenceTree.hpp"

using namespace std;

namespace {

SequenceTree makeThreeLeafTree() {
	istringstream newick("(Alpha:1,Beta:1,Gamma:1);");
	return SequenceTree(newick);
}

}  // namespace

// Single-line sequences (no interleaving needed) - the baseline case.
static void test_single_line_sequences() {
	SequenceTree tree = makeThreeLeafTree();

	istringstream phylip(
	    " 3 8\n"
	    "Alpha     ACGTACGT\n"
	    "Beta      TTTTAAAA\n"
	    "Gamma     GGGGCCCC\n");
	tree.mapSequencesOntoTree(phylip);

	SequenceTree::NodeVector nodes;
	tree.addNodesInPrefixOrder(nodes);

	string alpha, beta, gamma;
	for (auto *n : nodes) {
		if (NAME(n) == "Alpha") { alpha = SEQ(n);
}
		if (NAME(n) == "Beta") { beta = SEQ(n);
}
		if (NAME(n) == "Gamma") { gamma = SEQ(n);
}
	}
	assert(alpha == "acgtacgt");
	assert(beta == "ttttaaaa");
	assert(gamma == "ggggcccc");
}

// Forces the interleaved-continuation loop (phylipReadInterleavedContinuation)
// to actually run, exercising the migrated shared implementation rather
// than just the opening-line read. Continuation lines are padded with
// 10 leading spaces - see the file header comment for why plain,
// unpadded ones (this project's usual convention) don't work here.
static void test_multiline_interleaved_sequences() {
	SequenceTree tree = makeThreeLeafTree();

	istringstream phylip(
	    " 3 16\n"
	    "Alpha     ACGTACGTAC\n"
	    "Beta      TTTTAAAATT\n"
	    "Gamma     GGGGCCCCGG\n"
	    "          ACGTAC\n"
	    "          TTAATT\n"
	    "          GGCCCC\n");
	tree.mapSequencesOntoTree(phylip);

	SequenceTree::NodeVector nodes;
	tree.addNodesInPrefixOrder(nodes);

	string alpha, beta, gamma;
	for (auto *n : nodes) {
		if (NAME(n) == "Alpha") { alpha = SEQ(n);
}
		if (NAME(n) == "Beta") { beta = SEQ(n);
}
		if (NAME(n) == "Gamma") { gamma = SEQ(n);
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
static void test_unmatched_name_does_not_corrupt_real_sequences() {
	SequenceTree tree = makeThreeLeafTree();

	istringstream phylip(
	    " 5 16\n"
	    "Alpha     ACGTACGTAC\n"
	    "Beta      TTTTAAAATT\n"
	    "Fake1     AAAAAAAAAAAAAAAA\n"
	    "Gamma     GGGGCCCCGG\n"
	    "Fake2     CCCCCCCCCCCCCCCC\n"
	    "          ACGTAC\n"
	    "          TTAATT\n"
	    "          AAAAAA\n"
	    "          GGCCCC\n"
	    "          CCCCCC\n");
	tree.mapSequencesOntoTree(phylip);

	SequenceTree::NodeVector nodes;
	tree.addNodesInPrefixOrder(nodes);

	string alpha, beta, gamma;
	for (auto *n : nodes) {
		if (NAME(n) == "Alpha") { alpha = SEQ(n);
}
		if (NAME(n) == "Beta") { beta = SEQ(n);
}
		if (NAME(n) == "Gamma") { gamma = SEQ(n);
}
	}
	assert(alpha == "acgtacgtacacgtac");
	assert(beta == "ttttaaaattttaatt");
	assert(gamma == "ggggccccggggcccc");
}

int main() {  // NOLINT(bugprone-exception-escape) - a test binary crashing on an unexpected exception is a fine, loud failure mode.
	test_single_line_sequences();
	test_multiline_interleaved_sequences();
	test_unmatched_name_does_not_corrupt_real_sequences();
	cout << "SequenceTree_PhylipReader_test: all tests passed" << endl;
	return 0;
}
