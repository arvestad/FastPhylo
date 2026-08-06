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

### Phase B - Design the shared implementation (done)

**Scope refinement, decided while designing, not assumed away**:
Phase A's four differences all concern the *continuation-lines* phase
specifically (stop condition, `fin.eof()`, name-field-skip) - none of
them touch how each reader parses a sequence's *opening* line/name,
which uses three genuinely different mechanisms per site (whitespace-
delimited `fin >> names[i]` in `SequenceTree.cpp` vs. a fixed 10-char
`getline` field in the other two) that Phase A never investigated with
the same rigor. Consolidating *that* too would mean re-doing Phase A's
investigation for a second, un-analyzed set of differences - out of
scope for this plan as written. **This phase's shared implementation
covers only the continuation-lines loop and its completion check** -
the actual documented duplication (`anySequenceComplete()`/
`readInterleavedLines()`/`readInterleavedContinuationLines()`), not
the full top-level functions. Each reader keeps its own opening-line
logic untouched.

New header: `include/fastphylo/core/phylip_interleaved_reader.hpp`,
holding two small templates:

- `phylipAnySequenceComplete(numSequences, seqlen, getLength)` -
  replaces the three near-identical `anySequenceComplete()`
  definitions (and, per Phase A, `SequenceTree.cpp`'s single-pointer
  `actualNodeString` check too, since the any-of-all form is both
  equivalent under normal operation and safer).
- `phylipReadInterleavedContinuation(fin, numSequences, seqlen,
  getLength, appendLine)` - the shared "keep going until any position
  is done" loop. **Uses the corrected form Phase A proved out - no
  `|| fin.eof()`** - so the `DNA_b128_StringsFromPHYLIP` bug found in
  Phase A cannot recur in any of the three migrated call sites.

Deliberately *not* shared further: the per-position "read one line,
validate, append" step stays a callback (`appendLine(fin, i)`) that
each call site implements with its own mechanics (including that
site's own retry-on-blank-line and end-of-file handling) - the
mechanics genuinely differ (char-by-char validated reads vs.
`getline`-then-process through three different, differently-
permissive append functions), matching the same "don't force an
unnatural unification" judgment already applied to
`distance_matrix_refactor_plan.md`'s `fillMatrix_Hamming`.

The name-field-skip difference (`SequenceTree.cpp` only) stays local
to that reader's own `appendLine` callback, not parameterized into the
shared templates at all - it's specific to one call site's mechanics,
not a flag the shared loop itself needs to know about.

### Phase C - Migrate each call site, one at a time (done)

**`DNA_b128_StringsFromPHYLIP()`** (`Sequences2DistanceMatrix.cpp`) -
migrated first, since it's the one with the confirmed bug. Verified
via full rebuild + `ctest` + `RunExamples.sh` (byte-identical,
examples 1/2 exercise it directly), plus re-running Phase A's exact
bug reproduction (a constructed multi-line interleaved PHYLIP file
with and without a trailing newline): both variants now produce
byte-identical output, matching what `Sequence::readSequences` already
did correctly - the bug is fixed.

**`Sequence::readSequences()`** (`Sequence.cpp`) - migrated second.
This reader never had the `fin.eof()` bug, so behavior is unchanged,
only the duplication is gone. Verified via full rebuild + `ctest` +
`RunExamples.sh` (byte-identical, examples 6/11-18 exercise it via
`fastprot -I phylip`).

**`SequenceTree::mapSequencesOntoTree(std::istream&)`** - migrated
last, and needed the standalone-test approach since it has no live
caller. Written as a real, committed, CMake-wired test
(`tests/SequenceTree_PhylipReader_test.cpp`, linked against the
`fastphylo` library target directly and registered via `add_test()`),
not a throwaway - this is now permanent regression coverage for a
function that had none before.

**Two more things found and confirmed only by finally being able to
exercise this function**, neither visible from reading the code alone:

- `readSequenceLine()`'s nucleotide-normalization (`nucleotide2char()`)
  lowercases every character - a first draft of the test asserted
  uppercase expected strings and failed immediately; not a bug, just
  undocumented behavior worth having in the test's own comments now.
- **This reader's continuation-line convention is genuinely
  incompatible with the plain, unpadded continuation lines the rest of
  this project's PHYLIP files use** - the same convention
  `Sequences2DistanceMatrix.cpp`/`Sequence.cpp` read correctly and
  `RunExamples.sh`'s real fixtures use. Confirmed by direct experiment:
  a plain continuation line (matching this project's actual file
  format) throws `Bad character` mid-read, because the peek-based
  name-field-skip logic (difference #4) treats any continuation line
  starting with non-whitespace as having a stray 10-character name
  field to discard, silently eating real sequence data instead. The
  *same* content, padded with 10 leading spaces (as if every
  continuation line repeated a blank name field), reads correctly.
  This means `SequenceTree::mapSequencesOntoTree(istream&)`, if it
  were ever fed a real multi-line PHYLIP file in this project's normal
  format, would throw - it only works today (to the extent "today"
  means "if it were ever called at all") against a differently-padded
  PHYLIP dialect than the rest of the project uses. Whether its one
  real caller (`Simulator.cpp`, via externally-generated Seq-Gen
  output, itself only reachable from unbuilt `simulated_phylogenies/`
  tools) ever actually produces that dialect is unknown and out of
  this plan's scope to chase down - flagged here as a genuine,
  newly-discovered defect for whoever next touches that code path, not
  fixed as part of this consolidation (the name-field-skip logic
  itself is untouched, ported verbatim per Phase B's scope decision).
  The committed test's multi-line fixtures all use the verified-working
  space-padded form for this reason, documented in the test file's own
  header comment.

Also verified, empirically, that the `garbage`-exclusion sentinel
(Phase B's design) is load-bearing, not just theoretically motivated:
temporarily removed it and re-ran the test suite, confirming the
`test_unmatched_name_does_not_corrupt_real_sequences` test fails
exactly as predicted (`"Sequence not of correct length: Gamma"`) -
then restored the correct implementation and re-verified all three
tests pass, plus a full `ctest`/`RunExamples.sh` re-run.

### Phase D - Final verification (done)

Full clean rebuild (fresh `build/`): 0 errors. `ctest`: 3/3, including
the new `SequenceTree_PhylipReader_test` (permanent regression
coverage for a function that had none before this plan).
`RunExamples.sh`: byte-identical. `RunCliChecks.sh`: 61/61.

All 4 phases of this plan are now complete. Net result: the three
independently-copy-pasted interleaved-PHYLIP readers now share one
implementation of the part that was actually duplicated (the
continuation-lines loop and its completion check), a second real, live
bug is fixed (`DNA_b128_StringsFromPHYLIP`'s spurious exception on any
file without a trailing newline - same pattern as
`distance_matrix_refactor_plan.md`'s `fillMatrixRow_JC` finding), and
`SequenceTree::mapSequencesOntoTree(istream&)` has real test coverage
for the first time, which is also how its own genuine incompatibility
with this project's normal PHYLIP continuation-line format was found -
flagged for whoever next touches that code path, not fixed here.

## Decisions

- Branch: `phylip-reader-consolidation`, off `modernize-cpp17`'s tip
  (after merging `distance_matrix_refactor_plan.md`'s branch in first,
  since that plan also touched `Sequences2DistanceMatrix.cpp` and two
  diverged branches editing the same file independently was worth
  avoiding).
- Not blocking the rest of `lint_plan.md` (Phases 5's remaining two
  categories, or Phases 6-7) - this is follow-up work, flagged at
  review time, not scheduled ahead of the lint pass finishing.
