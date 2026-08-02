# Plan: consolidate the three interleaved-PHYLIP-sequence readers

## Goal

Three independent implementations of "read a PHYLIP file whose
sequences may be split across multiple interleaved lines" exist in the
codebase today:

- `SequenceTree::mapSequencesOntoTree(std::istream&)` in
  `src/fastphylo/core/SequenceTree.cpp`, backed by local helpers
  `readSequenceLine()` / `readInterleavedContinuationLines()`.
- `DNA_b128_StringsFromPHYLIP()` in
  `src/fastphylo/dna/Sequences2DistanceMatrix.cpp`, backed by local
  helpers `anySequenceComplete()` / `readInterleavedLines()`.
- `Sequence::readSequences(std::vector<Sequence>&, std::istream&)` in
  `src/fastphylo/core/Sequence.cpp`, backed by its own (same-named,
  separately-defined) `anySequenceComplete()` / `readInterleavedLines()`.

All three were touched during `lint_plan.md`'s Phase 5
(`readability-function-cognitive-complexity` - see
`lint_phase5_refactors.md`), each getting its own in-file extraction to
clear the complexity threshold. That pass deliberately stopped short of
merging them into one shared implementation - the three copies differ
in real ways (see below), and untangling that safely is its own piece
of work, not a lint-metric fix. This plan is that work.

## The three copies are variations on a theme, not identical twins

This is the part worth being explicit about, since it's exactly what
makes a naive merge unsafe:

1. **Different data structures.** `SequenceTree.cpp`'s version writes
   into `SequenceTree::Node`s, resolved through a name→node hash map
   built from the tree, with an explicit "garbage" bucket for names
   that don't match any node. The other two write positionally into a
   `std::vector<DNA_b128_String>` / `std::vector<Sequence>` with no
   external name resolution - the name read from the file *is* the
   record.

2. **Different stop conditions.** `SequenceTree.cpp` tracks one
   specific pointer - `actualNodeString`, the last node actually
   matched to a tree node during the first pass - and loops `while
   (actualNodeString->length() < seqlen)`. The other two instead scan
   *every* sequence each time via `anySequenceComplete()`, stopping as
   soon as **any one** of them reaches `seqlen` - a different
   mechanism, preserved exactly during Phase 5 rather than
   "corrected" to whichever reads more intuitively, since changing
   stop-condition semantics is not what a complexity-reduction pass
   should be doing.

3. **Even the two `anySequenceComplete()`-based copies differ from
   each other**: `Sequences2DistanceMatrix.cpp`'s outer loop is `while
   (whileTrue || fin.eof())`; `Sequence.cpp`'s is `while (whileTrue)` -
   no `|| fin.eof()`. Also preserved as-is during Phase 5.

4. **Only `SequenceTree.cpp`'s continuation-line loop skips a leading
   name field** (`peek()`, then skip 10 chars if the peeked character
   isn't whitespace). The other two read and append every continuation
   line unconditionally, with no such handling.

A safe consolidation has to either **prove** each of these differences
is a no-op restated differently (some might be - the `bootstrapSequences`
`stride`-chunking finding elsewhere in this same Phase 5 pass turned
out to be exactly that, proven via reasoning about `rand()` call order
before being unified across three files), or **parameterize** the
shared implementation over whichever ones turn out to be real.
Assuming either answer without checking is how this kind of merge
introduces a silent regression in a format-parsing path that has no
dedicated test coverage today.

## Why this is its own plan, not a lint-pass line item

Same shape of reasoning as `distance_matrix_refactor_plan.md`'s
`fillMatrix*` consolidation (also still not started): the individual
complexity fix is done and safe; the deeper duplication is real,
correctness-sensitive, and worth a dedicated look rather than a
drive-by decision made while chasing an unrelated metric. Specifically
here:

- The three call sites operate on genuinely different data structures,
  so "shared implementation" means real abstraction (a small
  interface/callback for "get sequence *i*'s current length" / "append
  to sequence *i*"), not just deleting two of the three copies.
- This is legacy, historically fragile parsing code - Phase 5 already
  turned up real subtlety while doing the mechanical extractions (see
  `lint_phase5_refactors.md`'s `SequenceTree::mapSequencesOntoTree`/
  `Sequences2DistanceMatrix.cpp`/`Sequence.cpp` sections), and none of
  the three has dedicated unit-test coverage - two are reachable only
  through `RunExamples.sh`'s end-to-end fixtures, one isn't reachable
  from any shipped binary at all (see below).

## What's in scope

- `src/fastphylo/core/SequenceTree.cpp`:
  `mapSequencesOntoTree(std::istream&)`, `readSequenceLine()`,
  `readInterleavedContinuationLines()`.
- `src/fastphylo/dna/Sequences2DistanceMatrix.cpp`:
  `DNA_b128_StringsFromPHYLIP()`, `anySequenceComplete()`,
  `readInterleavedLines()`.
- `src/fastphylo/core/Sequence.cpp`: `Sequence::readSequences()`, its
  own `anySequenceComplete()` / `readInterleavedLines()`.

**Reachability, for verification planning** (established during Phase
5, not re-derived here):

- `DNA_b128_StringsFromPHYLIP()` is called by `fastdist`'s
  `PhylipMaInputStream` - `RunExamples.sh`'s `ex1`/`ex2` exercise it.
- `Sequence::readSequences()` is called by the shared
  `io/PhylipMaInputStream.cpp`, used by `fastprot -I phylip` -
  `RunExamples.sh`'s `ex6`/`ex11`-`ex18` exercise it.
- `SequenceTree::mapSequencesOntoTree(std::istream&)` has **no live
  caller** in `fastdist`/`fnj`/`fastprot` - its only real caller is
  `src/c++/Simulator.cpp`, itself only reachable from
  `src/c++/simulated_phylogenies/` (unbuilt, not in any
  `CMakeLists.txt`). Verified during Phase 5 via a standalone test
  program compiled directly against `libfastphylo.a` (not committed) -
  the same approach will be needed here.

## Phases

### Phase A - Prove or disprove each of the four differences

For each difference listed above, determine definitively whether it's
a real behavioral distinction or a no-op restated differently, the
same rigor used for the `bootstrapSequences` stride-chunking finding
elsewhere in this Phase 5 pass (reasoned proof + empirical byte-for-byte
comparison, not just one or the other). Write up the finding for each
of the four before touching any code - this phase produces the actual
design constraints Phase B needs, rather than discovering them
mid-refactor.

### Phase B - Design the shared implementation

Once Phase A settles which differences are real, design one shared
implementation (an algorithm plus a small data-access
abstraction/callback covering "sequence *i*'s current length" / "the
buffer to append line *i* into"), parameterized over whichever of the
stop-condition and name-field-skip differences Phase A found to be
genuine. If a difference turns out to be a no-op, don't parameterize
over it - collapse it, same as the `stride`-chunking unifications
already done in this file family.

### Phase C - Migrate each call site, one at a time

Rebuild + `ctest` + `RunExamples.sh` after each migration for the two
reachable copies (`DNA_b128_StringsFromPHYLIP()`,
`Sequence::readSequences()`). `SequenceTree::mapSequencesOntoTree
(std::istream&)` needs the standalone-test-against-`libfastphylo.a`
approach from Phase 5 (write it as a committed test this time, given
it's the only way to verify this code path at all, rather than a
throwaway).

### Phase D - Final verification

Full clean rebuild + `ctest` + `RunExamples.sh`, plus a rerun of
whatever standalone test(s) Phase C produced for the unreachable copy.

## Decisions

- Branch: TBD - continue on `lint-cleanup`, or its own branch, given
  this is structural rather than lint-mechanical (same open question
  `distance_matrix_refactor_plan.md` has for its own scope).
- Not blocking the rest of `lint_plan.md` (Phases 5's remaining two
  categories, or Phases 6-7) - this is follow-up work, flagged at
  review time, not scheduled ahead of the lint pass finishing.
