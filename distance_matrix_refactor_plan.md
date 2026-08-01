# Plan: DistanceMatrix/FloatDistanceMatrix/DistanceRow index-type refactor

## Goal

`DistanceMatrix`, `FloatDistanceMatrix`, and `DistanceRow` (all in
`include/fastphylo/core/`, sharing a common macro-templated design -
`DISTANCEMATRIX`/`DM_TEMPLATE` etc.) index their data with `int`:

```cpp
DistanceType getDistance(int i, int j) const;
void setDistance(int i, int j, DistanceType d);
void setIdentifier(int i, Identifier id);
Identifier& getIdentifier(int i);
```

Nearly every caller today builds its indices from `.size()` or a loop
counter, both `size_t` by modern convention. That mismatch is the
direct cause of `lint_plan.md`'s Phase 2 finding 300+
`bugprone-narrowing-conversions` warnings, concentrated at these
exact call sites, across 15 files in the active build (and 8 more in
programs that are permanently out of scope for this engagement but
still compile against the same headers - see "What's in scope" below).

This isn't 300 independent minor issues - it's one API that predates
the rest of the codebase's shift to `size_t`-based indexing, echoing
out to every call site. `lint_plan.md`'s Phase 2 is explicitly
deferring to this plan rather than mechanically casting at each of
those 300+ sites, which would just re-endorse the mismatch instead of
fixing it.

## DistanceMatrix, FloatDistanceMatrix, and DistanceRow are copy-paste siblings, not a hierarchy

Worth being explicit about, since it changes what "fix the type" means
here: these are three **independently-written, structurally near-
identical class templates**, not one design reused three times.

- `DistanceMatrix` (Isaac Elias, 2006, the original) - full
  **symmetric** matrix, stores only the upper triangle (`i<=j ? D[i][j]
  : D[j][i]`).
- `FloatDistanceMatrix` (Mehmood Alam Khan, added Dec 2011 - five years
  later) - full **asymmetric** matrix, stores every cell. Its own
  header comment literally changed from "A symetric distance matrix"
  to "An A-symetric distance matrix" - used for the memory-efficient
  row-streaming output path.
- `DistanceRow` (Richard Schobesberger, undated) - a single row, same
  shape again, also for the streaming path.

All three independently define the *same-named* macros (`DM_TEMPLATE`,
`DISTVEC`, ...) in their own headers and share the identical `int i,
int j` accessor shape - the signature of copy-paste-derived siblings,
not shared code. This means the `int`/`size_t` mismatch and any `-1`
sentinel pattern are duplicated three times across three
independently-maintained classes, not one place - Phase A's inventory
must cover all three separately rather than assuming a fix to one
carries over.

## Why this needs to be its own plan, not a lint-pass line item

Tracing the narrowing warnings turned up something that makes a blind
"just widen the parameters to `size_t`" change unsafe to do
mechanically: at least one real, working piece of logic depends on
these fields being unsigned and wrapping around `-1`.

`Sequences2DistanceMatrix.cpp`'s `fillMatrixRow_{K2P,Hamming,JC,TN93}`
(four near-identical ~90-line functions, one per distance model - see
"The Sequences2DistanceMatrix.cpp duplication" below) each contain:

```cpp
} else {
    row = -1;                                    // row is size_t
}
for (size_t j = row + 1; j < numSequences; j++) { // -1 wraps to SIZE_MAX, +1 wraps back to 0
```

Confirmed (not assumed) that `row` is never passed to `getDistance`/
`setDistance` while holding that wrapped value in any of the four
copies - so this specific instance isn't a live bug - but it's exactly
the kind of fragile, non-obvious code an automated or careless
widening pass could break without anyone noticing (e.g. someone
"clarifying" `row = -1` to `row = 0` would silently skip index 0).

Two more instances of the same `field = -1` pattern on unsigned
members exist in the library itself, not just in caller code:

- `DistanceRow_impl.hpp:38`, in the `DistanceRow(istream&)`
  constructor - directly under a `//Doesn't function here :)` comment
  from the original author, i.e. already self-flagged as broken/dead,
  not a working trick. `columns = -1;` is followed immediately by
  `objInitFromStream(in)`, which presumably repopulates `columns`
  properly, making the assignment likely vestigial - needs
  confirming, not assuming, same as everything else here.
- `FloatDistanceMatrix_impl.hpp:40-41`, in the equivalent
  `FloatDistanceMatrix(istream&)` constructor - `rows = -1; columns =
  -1;` before `objInitFromStream(in)`. Same shape, no "doesn't
  function" disclaimer this time; needs its own check.

`fnj/XmlInputStream.cpp`'s `l.row_nr = -1` is a different, safe case -
`row_nr` is declared `int` (`XmlInputStream.hpp:19`), so no unsigned
wraparound is involved there. Noted so it isn't re-investigated later
as if it were the same risk.

**The point**: this codebase already relies on signed/unsigned
interplay as a working feature in at least one place, and has at
least two more `-1`-on-unsigned assignments of unconfirmed status.
Widening the API's parameter types without inventorying every such
site first would either (a) silently break the working trick, or (b)
be silently no-safer than before if a genuinely-live sentinel gets
missed. This plan's Phase A is exactly that inventory.

## The Sequences2DistanceMatrix.cpp duplication

Independent of the type question, `fillMatrix_{Hamming,JC,K2P,TN93}`
and `fillMatrixRow_{Hamming,JC,K2P,TN93}` (8 functions total, ~80-100
lines each) are near-identical: same skeleton (iterate over sequence
pairs, resolve ambiguities, fill the matrix/row), differing only in
which per-pair distance formula gets called in the middle. The `-1`
wraparound trick above is duplicated across all four `fillMatrixRow_*`
copies as a direct consequence of this duplication - fixing the
duplication fixes that four-times-over.

Proposed shape: one shared implementation per family (full-matrix and
row-streaming), each parameterized by the per-pair distance callable
(a function pointer or small `std::function`/template parameter,
matching what `compute_K2P`/`compute_Hamming`/`compute_JC`/
`compute_Tamura_Nei` already look like), collapsing 8 functions to 2.
The `row = -1` trick gets replaced with something a reader doesn't
have to reverse-engineer - e.g. an explicit `std::optional<size_t>
startRow` or a named `bool fromBeginning` parameter - as a natural
side effect of writing the shared version once instead of four times.

This refactor is *why* the type fix and the duplication cleanup belong
in the same plan: doing the type fix first without touching the
duplication means fixing the same wraparound trick four times over;
doing the duplication cleanup first makes the type fix trivial to
apply once, in the new shared function, instead of a call-site
still to be reconciled with the same design decision in four places.

## What's in scope

`include/fastphylo/core/{DistanceMatrix,FloatDistanceMatrix,
DistanceRow}.hpp` (+ their `_impl.hpp` bodies) and every caller in the
active build: `apps/{fastdist,fnj,fastprot}`, `src/fastphylo/{core,
dna,protein}`. Confirmed 526 total call-site occurrences of
`getDistance`/`setDistance`/`getIdentifier`/`setIdentifier` across 23
files; roughly 15 of those files are in the active build.

`protein/Matrix.hpp` is a different class, already indexed with
`std::size_t` - confirmed via grep, not affected by this issue, not in
scope here (it's the subject of a separate, already-noted, unrelated
future investigation - fixed-size/Eigen-backed representation - see
memory).

**Explicitly out of scope but affected by the API change**:
`fastprot_mpi/`, `iterative_tree_merger/`, `sequence_based_tree_
reconstruction/`, `buildtree.cpp` all call these same accessors but
are permanently excluded from this engagement's build/verification
(not built by default, `RunExamples.sh` doesn't exercise them). An
`int`->`size_t` signature change would compile-break them silently
with no CI in this environment to catch it. Phase B must still update
their call sites for consistency (a `git grep` + mechanical fix, same
motions as everywhere else) even though they can't be verified by
building - flag this explicitly when doing that phase rather than
skipping it and leaving them broken.

## Phases

### Phase A - Inventory, no code changes

Full enumeration of every `int`-indexed accessor across the three
classes and every call site (526 occurrences, 23 files). For each of
the three known `field = -1`-on-unsigned sites (`Sequences2Distance
Matrix.cpp`'s four `fillMatrixRow_*` copies, `DistanceRow`'s stream
constructor, `FloatDistanceMatrix`'s stream constructor), confirm by
tracing (not assuming) whether the wrapped value ever reaches an
accessor call, and classify each as: working-and-load-bearing (must be
preserved or deliberately replaced with an explicit named idiom),
dead/vestigial (can be deleted), or already-safe (int-typed, no
wraparound). Produces a concrete site-by-site list the rest of the
plan executes against, rather than discovering issues mid-refactor.

### Phase B - Widen the three classes' index parameters to `size_t`

`getDistance`/`setDistance`/`getIdentifier`/`setIdentifier` (and any
other `int`-indexed method Phase A's inventory turns up) across
`DistanceMatrix`/`FloatDistanceMatrix`/`DistanceRow`. Update every
caller in the active build, verifying via full rebuild + `ctest` +
`RunExamples.sh` after each file or logical cluster of files - same
discipline as every other phase of this engagement. Update the
out-of-scope callers (`fastprot_mpi/` etc.) too, on a best-effort,
unable-to-fully-verify basis, flagged the same way Phase C of
`project_layout_plan.md` flagged its own unverifiable `fastprot_mpi`
fix.

### Phase C - Consolidate `Sequences2DistanceMatrix.cpp`'s fillMatrix*/fillMatrixRow* families

Collapse the 4+4 near-identical functions into one parameterized
full-matrix implementation and one parameterized row-streaming
implementation, per "The Sequences2DistanceMatrix.cpp duplication"
above. Replace the `row = -1` wraparound trick with an explicit,
named idiom in the new shared function. This is real design work, not
mechanical - budget real time for getting the per-model distance
callable's signature right and verifying every one of the 4 (or 8)
original call sites still produces byte-identical output for its
model.

### Phase D - Remaining narrowing-conversions call sites

Everything Phase B doesn't touch: `compute_Tamura_Nei`/
`compute_Tamura_Nei_fixratio`'s `int strlen`/`numAs`/etc, which
already (per `lint_plan.md`'s Phase 2 investigation) must stay `int`
because of `strlen - sd.deletedPositions`'s signed subtraction - these
get per-call-site explicit casts, not a type change. Sweep for
whatever else remains after Phases B/C land, since some of the 352
original findings could turn out to be unrelated to the DistanceMatrix
family entirely.

### Phase E - Final verification

Full clean rebuild + `ctest` + `RunExamples.sh` + a benchmark re-run
(this touches `fillMatrix*`, which is on the hot path measured by
`bench_primitives`/`plan.md`'s speed2026a work).

## Decisions

- This plan supersedes `lint_plan.md`'s Phase 2 for the
  `bugprone-narrowing-conversions` category specifically - `lint_plan.md`
  will note the deferral rather than duplicate this plan's content.
- Do the duplication cleanup (Phase C) as part of this plan rather than
  filing it separately, since the type fix and the duplication are
  entangled (per "Why this needs to be its own plan" above) - fixing
  one without the other means redoing work.
- Branch: TBD - continue on `lint-cleanup`, or a dedicated branch.
  Given the scope (structural, not lint-mechanical), a dedicated
  branch off `lint-cleanup` or off `modernize-cpp17`'s tip is worth
  deciding before Phase A starts.
