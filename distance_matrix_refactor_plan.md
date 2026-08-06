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

### Phase A - Inventory, no code changes (done)

Full enumeration confirmed (not just the original estimate): **508
call-site occurrences across 28 files** (close to the original 526/23
estimate; the small drift is just codebase churn since that number was
written, not a sign the original count was wrong). Of those 28, the
in-scope set (active build) is 17 files; the rest are the already-known
out-of-scope callers.

**A fourth `field = -1`-on-unsigned site, missed by the plan's original
inventory**: `DistanceMatrix_impl.hpp:59`, in `DistanceMatrix(std::istream&)`
- `size = -1;` immediately followed by `objInitFromStream(in)`, which
unconditionally overwrites `size` before any accessor call (`if
(newSize != size) { size = newSize; ...}` - `newSize` is never
`SIZE_MAX` in practice, so this branch always fires). Same
"harmless, immediately overwritten" shape as the other two documented
sites below - traced, not assumed.

Classification of all four `field = -1`-on-unsigned sites, each traced
individually:
- `DistanceMatrix_impl.hpp:59` (`size = -1`) - **already-safe**:
  overwritten before any accessor call, confirmed above.
- `FloatDistanceMatrix_impl.hpp:40-41` (`rows = -1; columns = -1;`) -
  **already-safe**: `objInitFromStream()`'s own override unconditionally
  reassigns both fields before any accessor call, same shape as above.
- `DistanceRow_impl.hpp:38` (`columns = -1`) - **dead, and actually
  broken**, confirmed rather than assumed: `DistanceRow` never
  overrides `objInitFromStream()` (grep confirms no such override in
  `DistanceRow.hpp`), so this constructor calls `Object`'s no-op
  base version - `columns` stays at the wrapped `SIZE_MAX` forever
  and `D` is never resized. Matches the `//Doesn't function here :)`
  comment directly above it. Also confirmed **unused**: zero call
  sites for `DistanceRow(istream&)` anywhere in the active or
  out-of-scope source tree. Left in place, not deleted, since deleting
  a public constructor is a design decision outside this refactor's
  stated scope (index types), not something to fold in opportunistically.
- `Sequences2DistanceMatrix.cpp`'s four `fillMatrixRow_*` `row = -1`
  sites - re-confirmed safe during Phase C below (never reach an
  accessor while wrapped), and replaced with an explicit named idiom
  there per Phase C's own plan.

**`removeLastRow()`/`addRow()` underflow risk, noted but out of
scope**: both do unguarded `--rows`/`size--` on already-`size_t`
members (pre-existing, not part of this refactor's `int`/`size_t`
parameter question). Every in-scope caller (`NeighborJoining.hpp`/
`.cpp`, `LeastSquaresFit.cpp`) only calls `removeLastRow()` inside a
loop that stops at 3 remaining nodes, so the underflow is never
triggered in practice - confirmed via grepping every caller, not
assumed safe.

**The four out-of-scope files with call sites
(`iterative_tree_merger/`, `sequence_based_tree_reconstruction/`,
`buildtree.cpp`) are not wired into any CMake target at all** -
confirmed via grep, not assumed. This corrects the plan's original
premise that Phase B needed to "update them for consistency even
though they can't be verified by building": since they're not merely
unverifiable but **entirely uncompiled**, and `int`->`size_t` is never
a hard compile error (only implicit widening), there is no actual
compile-break risk to guard against by touching them. Phase B below
leaves them untouched rather than making blind, unverifiable edits to
dead code for stylistic consistency alone.

### Phase B - Widen the three classes' index parameters to `size_t` (done)

`getDistance`/`setDistance`/`getIdentifier`/`setIdentifier` widened
from `int`/mixed to `size_t` across all three classes
(`DistanceMatrix`, `FloatDistanceMatrix`, `DistanceRow`) - header-only
change, no `_impl.hpp` bodies needed touching since the accessor
bodies themselves already indexed `std::vector`s with the parameter
directly.

**Zero call-site changes were needed anywhere in the active build** -
confirmed by rebuilding immediately after the header change: 0
compile errors. Traced why beforehand (not just gotten lucky): every
in-scope caller already uses `size_t`-typed loop variables/`.size()`
results (verified directly in `NeighborJoining.hpp`, the single
largest caller at 136 call sites, all already `size_t` end to end) -
the `int` parameter was the one thing not already following the
codebase's own established convention, exactly as the plan's "Goal"
section framed it: one API predating the rest of the codebase's shift
to `size_t`, not 300+ independent problems.

Measured via a `clang-tidy --checks='-*,bugprone-narrowing-conversions'`
before/after diff across the same 57 in-scope files: **285 -> 170
findings, 115 resolved by this header-only change alone**, zero source
edits elsewhere. Verified: full rebuild (0 errors) + `ctest` (2/2) +
`RunExamples.sh` (byte-identical) after the change.

### Phase C - Consolidate `Sequences2DistanceMatrix.cpp`'s fillMatrix*/fillMatrixRow* families (done)

**Scope refinement, decided while implementing, not assumed away**:
the plan's original framing ("one parameterized full-matrix
implementation and one parameterized row-streaming implementation")
turned out to be exactly right for the row-streaming family but not
quite for the full-matrix family, once actually reading all 8
functions side by side rather than by their summary shape:

- `fillMatrixRowGeneric<RawDistance>()` - **all four**
  `fillMatrixRow_{Hamming,JC,K2P,TN93}` unify into this one template.
  All four originals had their ambiguity-resolution second pass
  entirely commented out already (dead code, deleted here rather than
  kept as an inert comment four times over) - once that's gone, even
  Hamming's row version (which in the full-matrix family below is
  genuinely different in shape) reduces to the exact same skeleton as
  the other three: compute a raw pairwise distance, convert it to a
  plain float, done.
- `fillMatrixML<RawDistance>()` - shared by `fillMatrix_JC`/`_K2P`/
  `_TN93` only. `fillMatrix_Hamming` stays its own function: its
  ambiguity handling is genuinely different in *shape*, not just
  formula - `compute_Hamming_distance()` returns a plain `float` (not
  an `ML_string_distance`), its closest-neighbor resolution calls the
  simpler `resolveAmbiguities()` (not
  `resolveAmbiguitiesUsingTransitionProbabilities()`), and its
  correction pass borrows `compute_JC()` as an ML estimate regardless
  of the model. Forcing Hamming into the same template as the other
  three would need enough extra parameterization (a
  "closest-neighbor-resolver" callable, a separate
  "correction-estimate" callable, a storage-shape switch between
  `RawDistance` alone and `pair<RawDistance,ML_string_distance>`) that
  the result would read *harder* than the duplication it replaced -
  judged not worth it. Net: 8 functions -> 3 pieces of source code
  (`fillMatrixML` instantiated 3x, `fillMatrix_Hamming` kept,
  `fillMatrixRowGeneric` instantiated 4x).

**A real, live, previously-uncaught bug found and fixed during this
consolidation**: `fillMatrixRow_JC`'s `mem_eff_flag` branch was
inverted relative to its three siblings (`if (mem_eff_flag)` instead
of `if (!mem_eff_flag)`) - and, since the condition governs *both*
branches via if/else, this meant JC's row-streaming path had the
*regular* (`--output-format=binary`) and *memory-efficient*
(`--memory-efficient`) behaviors swapped, not just one of them wrong.
Verified live, not just read: `fastdist -D JC -e -O phylip` on an
8-sequence test set printed zeros for the entire lower triangle
instead of the real distances (comparing against the full-matrix
reference, the correct behavior for `-e`), while the same command
with `-D K2P` produced output byte-identical to its own full-matrix
reference. Neither `RunExamples.sh` nor `RunCliChecks.sh` exercises
`fastdist`'s `--memory-efficient` flag at all, so this had zero
regression coverage before now. `fillMatrixRowGeneric()` uses the
three-out-of-four (correct) direction uniformly, so this specific bug
class cannot recur in any of the four models going forward.

Verified the fix didn't regress the three already-correct models nor
change unrelated behavior, via a genuine before/after comparison (a
`git worktree` at the pre-Phase-C commit, not just re-reading the new
code): `K2P`/`Hamming`/`TN93`'s binary-row-streaming output is
byte-identical before and after, for both `--output-format=binary`
and `--memory-efficient`. `JC`'s binary output changed by exactly one
bit (`-0.0` -> `0.0`, IEEE-754 negative vs. positive zero, numerically
and observably identical in every context) - traced to the fix
correctly skipping a now-redundant self-distance recomputation for
row 0 that the old buggy code's inverted branch used to redo. All
four models' full-matrix output (`fillMatrix()`, untouched by this
bug) is unchanged.

The `row = -1` wraparound trick (relying on unsigned integer overflow
to turn `row+1` into `0`) is gone from every model, replaced by an
explicit `startCol` variable in `fillMatrixRowGeneric()` - readable
without knowing wraparound is happening on purpose, and structurally
unable to invert silently the way JC's copy did. The
closest-neighbor/`closestDist` tracking in all four original row
functions was also dead (written, never read outside the same
commented-out block) - removed rather than carried into the new
shared code.

Verified: full clean rebuild (0 errors) + `ctest` (2/2) +
`RunExamples.sh` (byte-identical) + `RunCliChecks.sh` (61/61) +
`bench_primitives` (healthy speedups, 7x-68x depending on primitive/
length, no regression - this consolidation touches the same hot path
speed2026a measures) + the worktree-based before/after binary
comparison described above. Measured via the same
`clang-tidy --checks='-*,bugprone-narrowing-conversions'` diff used in
Phase B: 170 -> 150 findings (Phase B: 285 -> 170; this phase: a
further 20 resolved by removing the duplication, on top of the dead
code and wraparound trick that are gone regardless of what
`clang-tidy` reports).

### Phase D - Remaining narrowing-conversions call sites (done)

**Scope confirmed by direct discussion before starting**: of the 150
findings remaining after Phase C, most turned out to be a genuinely
different problem - 52+ were `double`->`float` precision loss (numeric
formula code, e.g. `TamuraNei.cpp`'s ts/tv-ratio math), plus stream/
`char` sentinel handling and libxml2 API calls, none of it about
index types. Explicitly scoped this phase to *only* the
DistanceMatrix-family-adjacent `size_t`->`int` sites, triaging each
remaining finding by hand (not a blind sweep) into three buckets:

- **Must-stay-`int`, cast at the call site** (the pattern the plan's
  own text already anticipated for `compute_Tamura_Nei`/
  `compute_Tamura_Nei_fixratio`'s `strlen`/`numAs`/`numCs`/`numGs`/
  `numTs` params - confirmed these must stay `int` due to
  `lint_plan.md`'s Phase 2 finding about signed subtraction inside
  them): applied `static_cast<int>(...)` at all 11 call sites this
  refactor's own Phase C code introduced in
  `Sequences2DistanceMatrix.cpp`. Found one more real instance of the
  same *class* of must-stay-signed constraint, not previously
  documented: `NeighborJoining.cpp`'s `computeNJTree()` counts `i`
  down through `0` with `i >= 0` - if `rows`/`columns` (fed by
  `dm.getRows()`/`getColumns()`, both `size_t` since Phase B) were
  ever widened instead of cast, the loop would never terminate (or
  read out of bounds) once `i` underflowed. Cast at the two call
  sites rather than widening, with a comment explaining why - the
  wrong fix here would have been a real, live bug, not just a lint
  regression.
- **A `size_t` value assigned into `str2int_hashmap`'s `int` value
  type** (`LeastSquaresFit.cpp`, `fnj/main.cpp`) - a project-wide
  utility type, not part of the DistanceMatrix family itself, so
  cast at the two call sites rather than changing the shared type.
- **Genuinely safe to widen** - the `DataOutputStream::printRow()`
  virtual hierarchy's `int row` parameter (base declaration +
  `BinaryDmOutputStream`/`PhylipDmOutputStream`/`XmlOutputStream`
  overrides, all four widened together to keep the override set
  matching - a mismatched override would silently stop overriding
  anything) and each implementation's `entriesPerRow` local variable,
  all simple upward-counting loop bounds fed directly by `size_t`
  values with no signed dependency. `PhylipDmOutputStream.cpp` had
  one call site already `static_cast<int>`'d then immediately
  `static_cast<size_t>`'d back for the loop condition - widening the
  variable outright removed both casts, not just silenced the
  warning. Same treatment for `applyFixFactor`/`applyFixFactorRow`
  (`DistanceMatrix.cpp`/`FloatDistanceMatrix.cpp`/`DistanceRow.cpp`) -
  their own free functions, directly part of this family, simple
  upward loops.

Explicitly left alone (confirmed unrelated by reading each site's
actual code, not assumed from the finding's category label alone):
`protein/Matrix.cpp`/`ExpectedDistance.cpp` (confirmed via grep - never
touch `DistanceMatrix`/`FloatDistanceMatrix`/`DistanceRow`/
`DistanceRow` at all, matches the plan's existing `protein/Matrix.hpp`
exclusion), `DNA_b128_String.cpp` (its own internal capacity/index
bookkeeping, unrelated), libxml2 API calls in both
`XmlInputStream.cpp` files (`xmlRelaxNGNewMemParserCtxt`'s own
signature, not project code), and `Sequences2DistanceMatrix.cpp`'s two
pre-existing PHYLIP-parsing/bootstrap-resampling findings (unrelated
to the fillMatrix family this plan is about).

Measured via the same `clang-tidy` diff, on a freshly regenerated
`compile_commands.json` after one run showed a spurious
"Found compiler error(s)" on `TamuraNei.cpp` from a stale
include-path cache (a tooling artifact, not a real code issue - the
actual `cmake --build` succeeded with 0 errors throughout, and
`TamuraNei.cpp` was never touched in this phase): **150 -> 111
findings**, all 39 resolved sites confirmed to be the
DistanceMatrix-adjacent `size_t`->`int` category; every remaining
finding re-confirmed to be the out-of-scope `double`->`float`/stream/
libxml2 categories, not a fix left undone.

### Phase E - Final verification (done)

Full clean rebuild (fresh `build/`, regenerated
`compile_commands.json`): 0 errors. `ctest`: 2/2. `RunExamples.sh`:
byte-identical. `RunCliChecks.sh`: 61/61. `bench_primitives`: healthy
speedups (4x-56x depending on primitive/length), no regression - this
phase's changes sit directly on the hot path `speed2026a`/
`bench_primitives` measures. Plus a DNA-distance sanity check across
all four models (`JC`/`K2P`/`TN93`/`HAMMING`) on a synthetic
300-seq/8000bp dataset, matching `lint_plan.md`'s own Phase 7
discipline - all fast (<0.1s) and consistent.

All 5 phases of this plan are now complete. The `int`/`size_t`
mismatch that motivated it is resolved everywhere it originated
(`DistanceMatrix`/`FloatDistanceMatrix`/`DistanceRow`'s own accessors)
and everywhere it was genuinely still visible in DistanceMatrix-
adjacent code; a real, previously-uncaught, live bug
(`fillMatrixRow_JC`'s inverted `--memory-efficient` behavior) was
found and fixed as a direct consequence of the Phase C consolidation
this plan called for from the start.

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
