#pragma once

#include <cstddef>
#include <istream>

// Shared by the three legacy interleaved-PHYLIP-sequence readers -
// SequenceTree::mapSequencesOntoTree(std::istream&) (SequenceTree.cpp),
// DNA_b128_StringsFromPHYLIP() (Sequences2DistanceMatrix.cpp), and
// Sequence::readSequences() (Sequence.cpp) - see
// phylip_reader_consolidation_plan.md for the investigation this is
// built from.
//
// Only the "keep reading interleaved continuation lines until any
// sequence reaches its target length" loop is shared here, not the
// per-line reading itself: the three readers' line-reading mechanics
// are genuinely different (SequenceTree.cpp reads and validates
// character-by-character via char2nucleotide(), throwing on a bad
// character; the other two `getline()` into a buffer first, then
// validate through their own, differently-permissive append
// functions), not just cosmetically different. Forcing that into a
// single shared function would need more parameterization than the
// duplication it replaced was worth (the same judgment call made for
// distance_matrix_refactor_plan.md's fillMatrix_Hamming).

// True once any of [0, numSequences) has reached seqlen, per
// getLength(i). Replaces the three near-identical
// anySequenceComplete() definitions this project had (one per reader)
// - and, for SequenceTree.cpp specifically, replaces its single-
// pointer actualNodeString check, which Phase A of the consolidation
// plan proved equivalent to this any-of-all form under normal
// operation (every position grows by exactly one line per pass, so
// checking one arbitrary matched position is equivalent to checking
// all of them) and additionally safer: actualNodeString starts
// nullptr and is only assigned if some name in the file matches a
// tree node, so the single-pointer form has a latent null-pointer
// dereference the any-of-all form cannot have.
template <class LengthFn> bool phylipAnySequenceComplete(int numSequences, std::size_t seqlen, LengthFn getLength)
{
    for (int i = 0; i < numSequences; i++)
    {
        if (getLength(i) == seqlen)
        {
            return true;
        }
    }
    return false;
}

// The interleaved-continuation-lines loop itself: repeatedly calls
// appendLine(fin, i) for every position until any one reaches seqlen.
// appendLine is responsible for that call site's own line-reading,
// validation, retry-on-blank-line, and end-of-file handling for a
// single position - the shared part is only "keep going until any
// position is done."
//
// Deliberately has no `|| fin.eof()` in its stop condition. Phase A
// of the consolidation plan proved that clause - present only in the
// pre-consolidation DNA_b128_StringsFromPHYLIP() - was a real, live
// bug: it made `fastdist -I phylip` throw a spurious exception on any
// file whose last line lacks a trailing newline (verified live with a
// constructed multi-line interleaved PHYLIP file; `fastprot -I
// phylip`, using the sibling reader without that clause, handled the
// identical file correctly). This is the corrected form, now shared
// by all three call sites - the bug cannot recur in any of them.
template <class LengthFn, class AppendLineFn>
void phylipReadInterleavedContinuation(std::istream &fin, int numSequences, std::size_t seqlen, LengthFn getLength,
                                       AppendLineFn appendLine)
{
    bool whileTrue = !phylipAnySequenceComplete(numSequences, seqlen, getLength);
    while (whileTrue)
    {
        for (int i = 0; i < numSequences; i++)
        {
            appendLine(fin, i);
            if (phylipAnySequenceComplete(numSequences, seqlen, getLength))
            {
                whileTrue = false;
            }
        }
    }
}
