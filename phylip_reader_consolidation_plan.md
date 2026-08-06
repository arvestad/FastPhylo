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

### Phase A - Prove or disprove each of the four differences (done)

Each of the four settled, by reasoning plus empirical test where the
reasoning alone wasn't conclusive - not assumed either way.

**1. Different data structures - confirmed real, not reconcilable.**
`SequenceTree.cpp` resolves each name through a hash map onto existing
tree nodes, with unmatched names redirected to a single shared
`garbage` string (cleared and reused per unmatched entry - confirmed
harmless, since nothing ever reads `garbage`'s content back).
`Sequences2DistanceMatrix.cpp`/`Sequence.cpp` write positionally into
a freshly-sized vector, no name resolution, every file entry kept.
This is inherent to what each caller needs (map onto an *existing*
tree vs. build a fresh list) - Phase B's shared implementation needs a
small data-access abstraction ("get me the buffer for position *i*"),
not a shared data structure.

**2. Different stop conditions - proven equivalent under normal
operation, but `SequenceTree.cpp`'s version has a latent crash the
other two don't.** `Sequences2DistanceMatrix.cpp`/`Sequence.cpp` check
`anySequenceComplete()` (any position reaching `seqlen`) each pass;
`SequenceTree.cpp` instead tracks one pointer,
`actualNodeString` (the last node matched during the first pass), and
checks only its length. Proof of equivalence: every position in the
vector-based readers grows by exactly one line per outer-loop pass
(confirmed by re-reading the inner loop - no path skips a position
without appending to it, aside from the empty-line retry which
re-attempts the *same* position rather than skipping ahead), so under
that uniform-growth invariant "any one reaches `seqlen`" implies "all
have" - checking one arbitrary matched position is equivalent to
checking all of them, and `SequenceTree.cpp`'s choice to track
specifically a *real* (non-`garbage`) position sidesteps `garbage`'s
own irregular growth rate (multiple unmatched names in the same file
would otherwise inflate a naive shared counter faster than the real
per-node rate). **But**: `actualNodeString` starts `nullptr` and is
only ever assigned inside the `if (iter != str2node.end())` branch -
if a PHYLIP file's names matched *zero* tree nodes, the continuation
loop's first condition check (`actualNodeString->length()`) would be a
null-pointer dereference. Not proven to be live (no reachable caller
exists to test end-to-end against - see "Reachability" above), but a
real latent defect independent of the consolidation question, and a
concrete argument for Phase B's shared implementation to use the
any-of-all-real-positions style uniformly rather than porting
`SequenceTree.cpp`'s one-pointer idiom forward - it's both equivalent
under normal operation *and* strictly safer.

**3. The `fin.eof()` difference between the two `anySequenceComplete()`
copies - confirmed a real, live, previously-uncaught bug, not a no-op.**
`Sequences2DistanceMatrix.cpp`'s outer loop is `while (whileTrue ||
fin.eof())`; `Sequence.cpp`'s is `while (whileTrue)`, no `fin.eof()`
clause. Verified live with a constructed multi-line interleaved PHYLIP
file (3 sequences, 250 chars each, forcing genuine multi-pass
continuation reading) in two variants, with and without a trailing
newline after the final data line:
- `fastdist -I phylip <file-without-trailing-newline>` (routes through
  `DNA_b128_StringsFromPHYLIP`) **throws immediately**: `"Sequence not
  of correct length: Alpha    length is 250"` - self-contradictory,
  since 250 *is* the correct `seqlen`. Root cause: the final content
  read (no trailing newline) sets `eofbit` on the *same* call that
  successfully completes the last sequence; `whileTrue` becomes
  `false` from `anySequenceComplete()`, but `fin.eof()` is now `true`,
  so `whileTrue || fin.eof()` is still `true` and the outer loop runs
  one more full pass over an already-exhausted stream, hitting the
  `myStr.empty() && fin.eof()` branch on its very first iteration.
- `fastdist -I phylip <file-with-trailing-newline>`: succeeds (the
  final content read stops at the newline delimiter without hitting
  real EOF, so `fin.eof()` is still `false` at the check).
- `fastprot -I phylip` (routes through `Sequence::readSequences`, no
  `fin.eof()` clause) succeeds identically on **both** file variants -
  byte-identical output regardless of trailing newline.

  This is a real, exploitable robustness gap in `fastdist`'s PHYLIP
  reader today: **any interleaved multi-line PHYLIP file whose last
  line lacks a trailing newline throws a spurious exception**, a very
  common real-world file-ending convention this project doesn't
  otherwise reject. Zero regression coverage exists for this
  (`RunExamples.sh`'s `seq.phylip` happens to end with a trailing
  newline). The shared implementation must use the
  `Sequence.cpp`/proven-correct form (no `fin.eof()` clause) -
  `DNA_b128_StringsFromPHYLIP`'s current behavior is the *wrong* one
  to preserve, not a variant to parameterize over. This is the second
  real bug this consolidation effort has turned up (after
  `distance_matrix_refactor_plan.md`'s `fillMatrixRow_JC` finding) -
  same pattern: duplicated code silently drifting into an actual
  defect in one copy, invisible until the copies are compared side by
  side.

**4. The name-field-skip, only in `SequenceTree.cpp`'s continuation
reader - confirmed real, and not safely droppable without more
information.** `SequenceTree.cpp`'s `readInterleavedContinuationLines()`
peeks at each continuation line's first character; if it's *not*
whitespace, it discards the next 10 characters before reading sequence
data, treating them as a stray name field a well-formed interleaved
file's continuation lines shouldn't have. Neither
`Sequences2DistanceMatrix.cpp` nor `Sequence.cpp` does this. Checked
whether it's actually load-bearing by reading how each of the three
readers' underlying append functions handle unexpected characters,
since that determines what "not doing this" actually costs:
`Sequence.cpp`'s `appendAllNonChars()` only strips literal `' '`
characters and blindly appends everything else, including
non-nucleotide letters, with **no validation at all** - a stray name
prefix would be silently spliced into the sequence data as if it were
real bases. `DNA_b128_String::append()` is pickier (a `switch` over
IUPAC nucleotide/ambiguity codes plus `' '`/`-`/`.`, falling through
to a loop-terminating `default` for anything else), so it would
*truncate* rather than corrupt - a different, but still real, failure
mode. Since `SequenceTree.cpp`'s function has no live caller to
determine whether real input ever actually exercises this path (see
"Reachability"), and the two other readers would each mishandle it
differently if it ever did occur, this isn't provably a no-op - it's a
genuine defensive feature specific to `SequenceTree.cpp`'s use case,
kept as an opt-in parameter in the shared implementation rather than
either dropped or forced onto the other two callers.

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
