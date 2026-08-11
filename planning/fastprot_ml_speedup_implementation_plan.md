# Plan: `fastprot` maximum-likelihood speedup implementation

## Relationship to the investigation

`planning/fastprot_ml_speedup_investigation_plan.md` was the
procedure - profile first, measure each of the four questions against
real data, don't reason from general library reputation. All five of
its phases are done. This document is the output: what to actually
build, in what order, and what not to bother with. Written 2026-08-11.

## What's already done

**Q3's finding shipped separately, ahead of this plan, because it
turned out to be a correctness bug, not a performance question**: PR
#2 (merged into `master`) replaced the finite-difference "Newton-
Raphson" solver with a safeguarded Newton search using analytic first
*and* second derivatives (ported from PHYLIP's `protdist.c`), fixing
a real bug where a large, model-dependent fraction of real-data pairs
(35.7%-94.3%) got unverified, sometimes badly wrong answers. Measured
side effect: 1.4-1.5x faster on well-behaved data, because true
Newton convergence (correct curvature) needs fewer iterations than
the noisy finite-difference approximation did, more than paying for
the one extra matrix multiply per iteration analytic curvature costs.

This is relevant to what follows in two ways:
1. It's already shipped, so the "current state" baseline for
   everything below is the *post-fix* code, not what Phase 1 profiled
   in early August 2026 (see "Re-baseline" below).
2. It changed the exact shape of the hot loop's per-iteration matrix
   work from two operations (one `MatrixExpm::at()`, one
   `Matrix::mult()`) to three (`MatrixExpm::at()` for `P(t)`, one
   `Matrix::mult()` for `P'(t)`, one more for `P''(t)`) - one more
   multiply than Q1's microbenchmark measured in isolation. Still the
   same operation shape (20x20 dense multiplies fed by a cached
   decomposition), so Q1's conclusion isn't invalidated, but the
   magnitude estimate below needs re-deriving, not just reused as-is.

## The one recommended change: Eigen for the ML hot path's matrix layer

**Do this. It's the only thing Phases 1-4 found worth doing.**

### What

Replace the `Matrix`/LAPACK-`dgemm_` implementation of the operations
`likelihood_slope_curv()` (`MaximumLikelihood.cpp`) performs every
Newton iteration - `MatrixExpm::at(t)`'s column-scale+multiply for
`P(t)`, and the two further `Matrix::mult()` calls for `P'(t)`/
`P''(t)` - with `Eigen::Matrix<double,20,20>` (fixed-size, stack-
allocated, no external BLAS call for a matrix this small).

**Scope precisely**: `MatrixExpm` and the specific `Matrix` operations
`likelihood_slope_curv()` uses. **Not** a wholesale `Matrix` class
replacement - `Matrix` has other callers (`ExpectedDistance.cpp`'s
`posterior_probability()`, `count_replacements()`'s tally-to-Matrix
conversion, general I/O) that aren't in this hot path and don't need
this. `Matrix`'s public API/ABI for those callers doesn't need to
change; only `MatrixExpm`'s internals and the handful of call sites
inside `likelihood_slope_curv()` do. Confirm this scope by grep before
starting - re-check `Matrix`'s caller list hasn't grown since Q1's
memory-based note about it (`project_matrix_implementation_investigation.md`).

### Why (Q1's evidence)

- Phase 1's profiling (pre-fix, but the same operation shape) found a
  single call site - the multiply inside `MatrixExpm::at()` - was
  45.1% of all profiled samples, going through Accelerate's general-
  purpose BLAS dispatch machinery for an 8,000-multiply-add 20x20
  problem it's not sized for.
- Q1's microbenchmark (`benchmarks/bench_matrix_backend.cpp`) measured
  the swap directly: Eigen `double` is 1.6-1.7x faster than the
  current path on the full column-scale+multiply operation (not just
  parity on the isolated multiply - the win is mostly Eigen fusing the
  scale-then-multiply into one allocation-free expression, where the
  current code allocates two `Matrix` temporaries per call). Eigen
  `float` is 3.2-3.6x faster. Both match the LAPACK-based result to
  far better than the 3-decimal-place tolerance (5e-16 and 1.7e-7
  respectively, against a 5e-4 bar).
- A naive hand-rolled loop is actively worse (~9x slower) - ruled out.
- MKL wasn't tested (a different BLAS vendor doesn't change the
  mechanism - dispatch overhead sized for large matrices - so it's not
  expected to help either) - not a candidate.

### Re-baseline before building (important, don't skip)

Q1's numbers were measured against the pre-solver-fix operation shape
(`MatrixExpm::at()` alone, and a separate raw `Matrix::mult()`), not
against the current three-operation `likelihood_slope_curv()`. A
straightforward extrapolation - summing Q1's per-operation numbers -
suggests the *whole* `likelihood_slope_curv()` call goes from roughly
2.6-2.8us (current) to roughly 2.0-2.3us with Eigen `double` (~1.2-
1.4x) or roughly 1.0-1.2us with Eigen `float` (~2.4-2.8x). **This is
an extrapolation, explicitly not a measurement** - the investigation's
own rule ("measurement, not general reasoning about which library
should be faster") applies here too. Before writing production code:
extend `bench_matrix_backend.cpp` (or write a fresh microbenchmark) to
time the actual current `likelihood_slope_curv()` end-to-end against
an Eigen-based reimplementation of the same three-operation function,
on the current (post-fix) code. Confirm the extrapolation roughly
holds before trusting it for a go/no-go decision on the `float` step
specifically (the `double` step is worth doing regardless, given the
naive-loop/MKL alternatives are already ruled out and `double` costs
nothing in accuracy).

**Done (2026-08-11)**, `benchmarks/bench_ml_hotpath.cpp`: measured
directly against the real, current `likelihood_slope_curv()`, all 5
models, real data (`examples/globin_family.fasta`). Confirms the
extrapolation's direction, a bit more conservative in magnitude:
**`double` ~1.20-1.24x, `float` ~2.06-2.13x**, consistent across all 5
models. Correctness: `double` agrees with the current LAPACK-based
result to 1e-13 to 1e-16 (machine precision); `float` to 1e-6 to 1e-8
in raw slope/curv units - three to five orders of magnitude inside
the solver's own `CONVERGENCE_TOL` (0.001), so `float`'s numerical
error is far too small to change which `t` the safeguarded-Newton
search converges to. Both steps confirmed worth building; proceeding
with `double` first per the sequencing below, `float` immediately
after given the margin here is comfortable enough not to need a
separate go/no-go pause.

### Steps

1. **Re-baseline microbenchmark** (above) - confirm the combined
   effect on the actual current hot-path shape, both `double` and
   `float`.
2. **Add Eigen as a real project dependency**: `vcpkg.json` (header-
   only, no new build complexity expected - already available via
   Homebrew in this dev environment, confirm `vcpkg`'s Eigen port
   works the same across the three CI platforms per
   `github_actions_release_builds_plan.md`'s cross-platform matrix).
3. **Verify `EigenSolver` vs. LAPACK `dgeev_` compatibility for real
   models, not just WAG** (Q1's microbenchmark only checked WAG) - run
   the same correctness check (max abs diff vs. the 5e-4 tolerance
   bar) for JTT/Dayhoff/MVR/LG too, and for a `-F`/`--model-file`
   custom matrix if a test fixture for that exists, since custom
   models are exactly the case Q1's design note flagged as needing
   the same verification, not assumed safe by analogy to WAG.
4. **Rewrite `MatrixExpm`'s internals and `likelihood_slope_curv()`'s
   three matrix operations** to use `Eigen::Matrix<double,20,20>`
   internally, converting to/from the existing `Matrix` type only at
   the function boundary (construction from `Q`, and the final
   `slope`/`curv` scalars) - `Matrix`'s public API elsewhere is
   unaffected. Byte-identical output is **not** expected (different
   eigensolver, different floating-point operation order) - use the
   3-decimal-place tolerance for verification, per Q1/Q3's established
   working tolerance, on `tests/MaximumLikelihood_test.cpp`'s existing
   real-data fixture (all 1500 pairs) plus the small fixture
   (`RunExamples.sh` example 18) and `RunCliChecks.sh`.
5. **Re-profile after landing the `double` step** (same method as
   Phase 1) to confirm the `dgemm_`/`libBLAS.dylib` cost is actually
   gone from the profile, the way `speed2026a`'s decomposition-caching
   fix was confirmed by its LAPACK-symbols-disappeared check.
6. **`float` as a separate, later, explicitly-decided step** - clear
   correctness margin per Q1, but decide with real re-baseline numbers
   in hand whether the further win over `double` alone is worth a
   second round of verification work, not bundled into step 4 by
   default.
7. **Update `INSTALL`/dependency docs** for the new Eigen dependency,
   matching how libxml2's optional-dependency documentation was
   handled in the `xml_optional_phylip_default` work.

### Explicitly not doing

- **MKL** - Q1's reasoning (dispatch overhead is vendor-agnostic at
  this matrix size) plus the real cost of provisioning it on all three
  CI platforms for an unlikely win.
- **Embedded eigendecomposition constants for the named models (Q2)**
  - measured at ≤0.3% of wall time even in the best case, smaller than
    run-to-run timing noise. Not worth the source verbosity or the
    maintenance burden of re-deriving constants if a model's published
    matrix is ever corrected.
- **Further precomputation sweep (Q4)** - done as part of writing this
  synthesis: re-read `calculate_ml_dists()`/`calculate_ed_dists()`
  (`ProtDistCalc.cpp`) fresh. Both already scope `Q`/`MatrixExpm`/
  `eq` construction to once per run, matching what Q2 already found
  and decided against optimizing further. No new candidates turned up.
- **A wholesale `Matrix` class rewrite** - `ExpectedDistance.cpp` and
  other non-hot-path callers don't need this; scope creep here adds
  risk (a class used in more places is more places to verify) for no
  measured benefit outside the ML hot path.

## `fastphylo-py` connection, revisited (per the investigation's own note to do so)

Phase 0 found `fastphylo-py` already does its own eigendecomposition
in Python (numpy/scipy) and only crosses into C++ for a thin, Brent's-
method-based branch-length optimizer taking pre-decomposed
eigenvectors/eigenvalues as plain arrays. This plan's Eigen work
doesn't change that picture much - `fastphylo-py` isn't calling
`fastprot`'s `Matrix`/`MatrixExpm` today regardless of what backend
they use internally, so an Eigen swap here doesn't by itself get
`fastphylo-py` anything. What *would* matter, if bindings are ever
pursued (still explicitly out of scope - "investigate and report
back," not a commitment): the decomposition-in/optimize-out boundary
`fastphylo-py` already assumes is the same shape this plan's
`MatrixExpm`/`likelihood_slope_curv()` split has. If `fastprot`'s
optimizer core (the safeguarded-Newton search from PR #2, now running
on Eigen types internally) were ever exposed as a small, self-
contained function taking a decomposition + observed counts and
returning a distance, `fastphylo-py` could in principle call *that*
instead of its own standalone `compute_protein_distances_cpp` - one
optimizer core shared by both, instead of two independent
implementations that could silently drift apart (as they already
partly have: `fastphylo-py`'s vendored core is stale, per Phase 0).
Not a recommendation to build this now - flagging it as the concrete
shape a future bindings effort should take, since it falls out
naturally from work this plan is already doing for other reasons.

## Sequencing

1. Re-baseline microbenchmark (cheap, answers the go/no-go on `float`
   before committing to it). **Done** - see above.
2. `double` Eigen swap (steps 2-5 above) - the actual work, verified
   change, real measured win. **Done** (`eigen-ml-matrix-backend`
   branch, off `master`): `MatrixExpm`'s internals and
   `likelihood_slope_curv()`'s three matrix operations now use
   `Eigen::Matrix<double,20,20>`; `MatrixExpm`'s public API unchanged
   (`Matrix` in/out via `at()`, used by the ED path); a new
   `at_eigen()`/`Q_eigen()` pair is the fast path
   `likelihood_slope_curv()` uses directly. `Eigen3` added to
   `vcpkg.json`/CI/`INSTALL` (`PUBLIC` link, since `Matrix.hpp` itself
   now needs Eigen's headers to declare `MatrixExpm`'s private state).
3. Re-profile, confirm. **Done** - macOS `sample`, 1000-sequence
   dataset: zero `DGEMM`/`DGEEV`/`libBLAS.dylib` samples (Phase 1's
   45.1% bottleneck is gone), replaced by Eigen's own GEMM kernels
   inlined into `fastprot` directly. **Real end-to-end wall-clock
   measurement** (`fastprot -D WAG -m`, interleaved reps, pre- vs.
   post-swap binaries, byte-distinct verified): **1.23x at 100
   sequences, 1.32x at 300, 1.35x at 600 and 1000** - increasing
   toward the larger, stated-hard-case sizes as expected (more of the
   run's wall time is inside the sped-up function as `N` grows),
   consistent with the microbenchmark's ~1.2x isolated-function
   prediction. Verified via full rebuild + `ctest` 5/5 on Clang
   Release and real GCC 14 (byte-identical `fastprot` output between
   compilers, all 5 models, real data), `-DWITH_LIBXML=OFF`,
   `RunExamples.sh` byte-identical on every fixture including the ML
   one, `RunCliChecks.sh` clean.
4. Decide on `float` with real numbers in hand. **Done** (2026-08-11),
   bundled with two further optimizations found while reviewing the
   code together (Lasse's questions directly led to items b and c
   below):
   a. **`float` for the per-iteration hot path.** `MatrixExpm` now
      caches a second, `float`-precision copy of the decomposition
      (`at_eigen_f()`/`Q_eigen_f()`) alongside the `double` one - the
      decomposition itself stays `double`-only (one-time cost, the
      more numerically sensitive step); only the cheap per-`t`
      evaluation runs in `float`.
   b. **`Q²` precomputed once per run** (`Q2_eigen_f()`), so
      `P'(t)=P(t)·Q` and `P''(t)=P(t)·Q²` are both computed directly
      from `P(t)` instead of `P''(t)` depending on `P'(t)` first - same
      multiply count, no dependency chain between the two.
   c. **`N` converted to a fixed `Eigen::Matrix<float,20,20>` once per
      pair**, before the Newton loop starts, instead of being re-read
      via `Matrix`'s bounds-checked `operator()` on every iteration
      (400 accessor calls/iteration, for a matrix that never changes
      across iterations - a real, previously-missed overhead. caught
      by Lasse asking "why is the hand-rolled Matrix class still used
      at all"). `likelihood_slope_curv()`'s reduction loop is now
      Eigen array expressions (`.array()`, `.cwiseMax`/`.max()`,
      `.square()`, `.sum()`) instead of a hand-written branchy loop -
      the branch (skip zero `N` cells) is gone; letting Eigen
      vectorize the whole 20x20 array uniformly was measured to not
      need it. The final reduction still accumulates in `double`
      (`.cast<double>()` before `.sum()`) - per-element arithmetic is
      where `float`'s speed comes from, summing 400 terms costs nothing
      extra in `double` and avoids adding the reduction's own rounding
      error on top of the already-comfortable margin below.

   **Considered and dropped**: exploiting `N`'s sparsity to compute
   only the `P'(t)`/`P''(t)` cells actually needed (skip full dense
   rows for amino acids that never appear as a "from" state).
   Skepticism from Lasse ("the non-zero elements can be just about
   anywhere") turned out to be right - measured directly
   (`benchmarks/measure_N_sparsity.cpp`) on both real
   (`examples/globin_family.fasta`) and synthetic data: `rows(N)`
   (distinct amino acids appearing as a "from" state) is essentially
   always all 20 (median 20 on both datasets, minimum 18 even on the
   smaller real one). Working through the FLOP count with that
   constraint gives at best ~1.3x even before accounting for a
   hand-rolled sparse loop losing to a vectorized dense kernel in
   practice - not worth the complexity. Kept the measurement script
   as a record of why, in case the same idea comes up again.

   **Measured** (interleaved reps, byte-distinct binaries verified,
   all 5 models, both 300 and 600-sequence datasets): **~2.0x-2.2x
   faster than the original (pre-Eigen, PR #2 baseline) implementation,
   ~1.5x-1.7x faster than the `double`-only Eigen step (item 2/3
   above)** - consistent across every model, holding from 300 to 600
   sequences. Correctness: `tests/MaximumLikelihood_test.cpp` (1500
   real pairs, all 5 models) passes unchanged. `RunExamples.sh`'s
   `ex18.out` shifted by ~1e-6 (regenerated, disclosed) - `float`'s own
   precision limit, as expected. New, disclosed property: Clang/GCC
   output is no longer guaranteed *byte*-identical with `float` in the
   hot path (checked: WAG showed single-last-digit differences of
   ~1e-6 between compilers on real data, JTT/DAY/MVR/LG were clean in
   this run but aren't guaranteed to always be) - ordinary cross-
   compiler floating-point codegen variance at `float`'s own precision
   edge (`exp()`/FMA differences), not a correctness issue, and still
   3+ orders of magnitude inside the solver's tolerance. This is a
   real, if narrow, departure from the `double` path's guaranteed
   bit-for-bit cross-compiler identity - worth knowing, not a defect.
   Verified via full rebuild + `ctest` 5/5 on Clang Release and real
   GCC 14, `-DWITH_LIBXML=OFF`, `RunExamples.sh`/`RunCliChecks.sh`
   clean, re-profiled (still dominated by Eigen's own GEMM kernel, now
   operating on `float` - no unexpected new hotspot).
5. Update dependency docs, changelog entry (`2.0.0-beta.4` or
   whatever's next, per `RELEASING.md`'s established prerelease
   convention). Dependency docs (`INSTALL`, CI workflows, `vcpkg.json`)
   done as part of step 2. Changelog entry and version bump: not done
   yet, pending PR review/merge.

## Side finding (2026-08-11): a real JTT data bug, found by benchmarking against PHYLIP

Lasse asked for a wall-clock comparison against PHYLIP's `protdist`
(built from the real 3.65 source - the same mirror used earlier this
engagement to study `makedists()`'s algorithm). Results, all JTT
model, same alignments, interleaved:

| N sequences | `protdist` | `fastprot` | speedup |
|---|---|---|---|
| 25 (real, `globin_family.fasta`) | 0.092s | 0.034s | 2.7x |
| 100 | 2.130s | 0.062s | 34.2x |
| 300 | 19.048s | 0.292s | 65.2x |
| 600 | 76.714s | 1.072s | 71.6x |

Real, but not solely a credit to this plan's work - `protdist` is
unoptimized reference C from 1993-2004, not a tuned implementation, so
a large gap is expected on its own. Worth having the number anyway,
since it's the standard tool people actually compare against.

**A spot-check of the actual distances (not just timing) found a real
bug**, not just a modeling-convention difference: JTT distances
between `fastprot` and `protdist` on the same pair differed by
~30-40%, non-uniformly across pairs (not a constant scale factor,
which would suggest a units/convention difference rather than a data
error). Traced to `ModelMatrix.cpp`'s `get_jtt()`: `Q(0,0)` was
hardcoded as `-1.24783051471057305` - **exactly 100x** the correct
value. Confirmed precisely: the sum of row 0's 19 off-diagonal entries
is `0.01247830514710573`, matching the stored (wrong) value divided
by 100 to 13 significant figures - a decimal-point transcription
error, not a modeling choice. This broke the fundamental "rows sum to
zero" invariant every valid continuous-time Markov rate matrix must
satisfy, corrupting `JTT`'s eigendecomposition (one spurious dominant
eigenvalue, -1.248, vs. everything else clustered around -0.001 to
-0.02) and therefore every JTT maximum-likelihood distance `fastprot`
has ever computed. `benchmarks/check_row_sums`-style sweep (see
`tests/ModelMatrix_test.cpp` below) confirmed this was an **isolated,
single-cell bug** - no other row of JTT, and no row of any other
model (`WAG`/`Dayhoff`/`MVR`/`LG`), had the same problem.

Fixed in both the main library (`src/fastphylo/protein/ModelMatrix.cpp`)
and `fastprot_mpi`'s mirrored, unbuildable-here copy (same reasoning
as the `ARVE` removal - a trivial, low-risk, mechanically-verifiable
data fix, worth keeping in sync even though the performance work isn't
mirrored there). After the fix, `fastprot`'s JTT distances agree with
`protdist`'s to within 0.41% (down from ~30-40%) - consistent with two
independent, correct implementations differing only in minor
precision/convergence details, not a real disagreement.

**Follow-up requested and done**: cross-validate all embedded model
matrices against IQ-TREE2's source (`model/modelprotein.cpp`, embedded
NEXUS-format lower-triangular exchangeability matrices + frequency
vectors, fetched and parsed - see `benchmarks/compare_to_iqtree2.py`
and `benchmarks/dump_model_matrices.cpp`). Built `Q_ref(i,j) =
R(i,j)*pi(j)` from IQ-TREE2's data, compared element-by-element
against `fastprot`'s matrices after solving for a single global scale
factor per model (empirical rate matrices are only defined up to an
arbitrary time-unit, so a per-model constant multiplier is expected
and fine - a per-*element* scale difference would not be, and is
exactly the kind of thing this check would catch). Result: **`WAG`,
`JTT` (post-fix), `Dayhoff`, and `LG` all agree with IQ-TREE2 to
0-0.003% (floating-point rounding of independently-transcribed
literals, not a real difference)**. `MVR` (Müller-Vingron) isn't in
IQ-TREE2's model catalog - no reference available there, but its rows
already passed the sum-to-zero check.

**Permanent regression coverage added**: `tests/ModelMatrix_test.cpp`
- checks every model's `Q` matrix rows sum to zero and equilibrium
frequencies sum to ~1, for all 5 models. This is the check that would
have caught the JTT bug immediately, with no external dependency -
`tests/MaximumLikelihood_test.cpp`'s existing check (does the solver's
returned answer satisfy its own convergence criterion) could not have
caught it, since a corrupted-but-still-valid-looking optimization
problem converges to *an* answer, just the wrong one; nothing checks
that answer against an external ground truth. No existing
`RunExamples.sh` fixture exercised JTT with `-m` at all, which is also
why the byte-identical-output checks never caught this.

## Follow-up (2026-08-11): closing the loop on "are we done"

Lasse pushed back on two points before treating this as finished -
both real gaps, not just process:

1. **The JTT speedup number was measured against the buggy matrix on
   the "old" side of the comparison** (`master` never got the data
   fix, only this branch did) - meaning part of the measured
   difference could have been the *data* changing, not just the code.
   Isolated properly (patched *only* the one-line JTT fix onto
   `master`, nothing else, so both sides solve the identical
   mathematical problem): **2.14x**, matching the originally-reported
   2.05x (measured against the buggy matrix) closely enough that the
   bug hadn't meaningfully distorted the headline number - but this
   should have been checked before reporting it, not after being
   asked.
2. **The bug also corrupted the non-ML (`ED`, no `-m`) path**, not
   just maximum likelihood - both call the same `get_jtt()`. Confirmed
   directly (same pair, `examples/globin_family.fasta`, JTT, no `-m`):
   `1.473649` (buggy, several entries pinned to a `0.001667` floor
   value - a clear sign of a badly broken optimization) vs. `2.047919`
   (fixed). Both now correct from the same single fix; noted here
   since it broadens the bug's real-world impact beyond what the
   original writeup described.

A fresh profile at this point (JTT, 1000 sequences) surfaced one more
thing worth fixing: `likelihood_calc()`'s early-exit check computes
`N.sum() - N.sum_diag()`, and `kimura_distance(N)` right below it
computes the *same* two aggregates over the *same* `N` again - 8.5% of
sampled time in a `Matrix::sum()`/`sum_diag()` call that only exists
because of this duplication. Lasse asked for three things in response:
fix the redundancy, reconsider whether `N`'s `Matrix` round-trip
conversion (`Matrix` in `ProtDistCalc.cpp` -> `Eigen` in
`likelihood_calc()`) is doing more work than it needs to, and actually
test raising `EIGEN_UNROLLING_LIMIT`.

**Redundancy + `N` conversion, done together** (they turned out to be
the same fix): `kimura_distance()`'s signature changed from taking `N`
to taking the two aggregates directly (`n_sum`, `n_sum_diag`) -
`likelihood_calc()` now computes them once (via Eigen's own
unchecked `.sum()`/`.trace()`, not `Matrix::sum()`/`sum_diag()`) and
reuses them for both the early-exit check and the starting value.
Separately, `N` itself no longer round-trips through the general
`Matrix` class at all for the ML path: `ProtDistCalc.cpp` gained
`tally_to_eigen()` (parallel to the existing `tally_to_matrix()`,
which the `ED` path still uses), converting `ProtSeqCode::
count_replacement_tally()`'s flat output directly into
`Eigen::Matrix<float,20,20>` - one conversion loop (no heap
allocation, no bounds-checked writes) instead of two (raw tally ->
`Matrix`, heap-allocated and bounds-checked, then `Matrix` -> `Eigen`,
bounds-checked reads, which is what `likelihood_calc()` used to do
internally on every call). `likelihood_calc()`/`likelihood_slope_curv()`
both now take `Eigen::Matrix<float,20,20>` directly - `Matrix` no
longer appears anywhere in `MaximumLikelihood.hpp` at all.

Measured (interleaved, byte-distinct binaries verified): **a further
1.16x-1.17x**, consistent across the two models checked (WAG, JTT).
Re-profiled after: the `Matrix::sum`/`sum_diag` bucket is gone from
the profile entirely (0 samples, down from 8.5%), confirming the fix
landed where expected.

**`EIGEN_UNROLLING_LIMIT` experiment** (`benchmarks/bench_unrolling.cpp`):
built as two separate standalone binaries (default limit 110 vs. a
raised 20000, comfortably above the ~8000 multiply-adds a 20x20x20
product costs) - two full binaries, not two code paths in one, since
mixing different values of this macro across translation units linked
into a single binary would be an ODR violation (the same template
instantiation would have two different bodies). **Result: no benefit.**
The isolated 20x20 multiply's median time was identical (0.333us) at
both settings - raising the macro didn't change the generated code for
this operation at all. Not adopted; kept the benchmark as a documented
record of a genuinely negative result, so this doesn't get re-tried
from scratch later.

Verified (all of the above): full rebuild + `ctest` 6/6 on Clang
Release and real GCC 14, `RunExamples.sh`/`RunCliChecks.sh` clean (no
fixture changed - this round was a pure refactor, output unchanged).

## Removed the internal "PAM" unit split (2026-08-11)

Lasse noticed the `100`/`0.01` split directly: `kimura_distance()`
computed in "PAM" units (`-100 * log(...)`, the historical convention
the published WAG/JTT/Dayhoff/MVR/LG matrices are calibrated in -
expected substitutions per 100 residues), and `calculate_ml_dists()`
(`ProtDistCalc.cpp`) rescaled the result by `0.01` at the very end to
get the substitutions-per-site units fastprot actually outputs (and
PHYLIP uses natively, with no internal PAM step at all). Preference:
drop the internal PAM representation entirely and work directly in
substitutions/site throughout, matching PHYLIP's own convention -
"an old idea that should go away."

The two scales aren't just a relabeling - `Q`'s eigenvalues are
calibrated for the *old* `t` (PAM units), so simply renaming `t`
without rescaling `Q` would silently break `e^{Qt}`'s meaning as a
transition-probability matrix. Rather than rescaling `Q` itself
(`ModelMatrix.cpp` - shared with the `ED`/non-`-m` path, which does
its *own*, separate `0.01 *` scaling and was deliberately left alone,
out of scope), `MatrixExpm`'s constructor gained an optional
`time_unit_scale` parameter (default 1.0, so every other caller -
`Matrix::expm()`/the `ED` path, tests using the default - is
unaffected): `m_Q = to_eigen20(Q) * time_unit_scale`, computed before
anything else, so the decomposition, its `float` copy, and `Q²` are
all automatically correct for the rescaled time unit with no further
special-casing. `calculate_ml_dists()` now constructs `MatrixExpm
Qdecomp(Q, PAM_TO_SUBSTITUTIONS_PER_SITE)` (`= 100.0`, named and
commented at its one call site) instead of the bare `MatrixExpm(Q)` -
every `t` flowing through `MaximumLikelihood.cpp` from that point on
is already in substitutions/site, with nothing left to rescale at the
output boundary (`ProtDistCalc.cpp`'s `0.01 * likelihood_calc(...)`
became a plain `likelihood_calc(...)`).

Every constant that depended on the old unit needed re-deriving, not
just dividing by 100 uniformly - `MAX_DISTANCE`/`MIN_DISTANCE` are
*distances* (500/1 -> 5.0/0.01, straightforward), but
`CONVERGENCE_TOL` is a threshold on the log-likelihood's *slope*
(1/distance-units): since `d(logL)/d(t/100) = 100 * d(logL)/dt`, the
equivalent threshold in the new units is **100x the old value**
(0.001 -> 0.1), not 100x smaller - the Newton step formula itself
(`t -= slope/curv`) turned out to need no special-casing at all, since
`curv` (a second derivative, 1/distance² units) rescales by 100² and
the ratio `slope/curv` naturally comes out in the right units either
way. `kimura_distance()`'s own `-100 * log(1 - adjusted)` became
`-log(1 - adjusted)` (the `t <= 0` fallback, previously the bare
literal `1`, now explicitly reads `MIN_DISTANCE` - was silently
assuming a value equal to a constant it never referenced).

One-line side benefit: the "too diverged" warning now prints the
*actual* clamped value (`distance=5`) instead of the old
internal-units number (`distance=500`, meaningless without knowing
about the `0.01` conversion that used to happen afterward) - verified
directly on a synthetic maximally-diverged pair.

Verified: full rebuild + `ctest` 6/6 on Clang Release and real GCC 14;
`RunExamples.sh` clean after regenerating `ex18.out` (shifted by
~1e-6, the same float-reordering-noise magnitude already characterized
for the `float` hot path - not a scaling bug); direct re-check against
PHYLIP's `protdist` on the same JTT pair gave the *exact* pre-refactor
value (`2.008679`), confirming this is a pure internal-representation
change with the real behavior unchanged.

## Fixed a second, independent scale bug: LG's rate matrix (2026-08-11)

While reviewing Dayhoff's matrix (its many zero entries are correct -
a known artifact of the small 1978 training dataset, and exactly why
JTT/WAG/LG were developed later; confirmed against IQ-TREE2's own
published Dayhoff matrix, same zeros), Lasse asked whether aligning
with IQ-TREE2's model data could sidestep any remaining scale
questions. That prompted computing each model's actual "mean
substitution rate" (`-sum(eq[i] * Q(i,i))`, the standard
IQ-TREE/RAxML/PAML convention for what a branch length of 1 means)
directly from IQ-TREE2's raw, unmodified `model/modelprotein.cpp`
data: WAG=95.28, JTT=100.64, Dayhoff=100.58, MVR has no IQ-TREE2
reference, but **LG=1.000001**.

This is a real, pre-existing bug, independent of the PAM-unit-removal
work above: the fixed `PAM_TO_SUBSTITUTIONS_PER_SITE = 100.0` constant
that `calculate_ml_dists()` used as `MatrixExpm`'s `time_unit_scale`
assumed every model's rate matrix was calibrated to the historical
"PAM" convention (mean rate ~100 - "100 expected substitutions per 100
residues" per unit time). That's true for WAG/JTT/Dayhoff/MVR, but LG
is published already normalized to a mean rate of 1 (the modern,
direct substitutions/site convention) - so scaling it by 100.0 dividing
gave `Qdecomp` calibrated 100x too coarse a time unit, producing ML
distances roughly two orders of magnitude too small for LG specifically
(caught via `benchmarks/analyze_convergence.cpp`: LG was the only model
hitting `MIN_DISTANCE` on 21-43% of real pairs, with far worse
convergence precision than the other four).

Fix: replaced the fixed constant with a new `mean_substitution_rate(Q,
eq)` function (`ModelMatrix.hpp`/`.cpp`) that computes a model's *own*
true mean rate from its data, so `calculate_ml_dists()` now passes
`1.0 / mean_substitution_rate(Q, eq)` as `time_unit_scale` - this
self-corrects for whatever convention a given model happens to be
published in (including any added later) rather than assuming one
fixed factor works for all of them. Also more numerically precise for
WAG (true rate ~95.3, not exactly 100) and MVR (true rate unknown but
not exactly 100 either) even though those two were never as badly
wrong as LG.

Re-ran `analyze_convergence.cpp` post-fix on both the globin family
(300 pairs) and the synthetic 300x300 benchmark (44,850 pairs), all 5
models: **100% converged, 0% hit either boundary, max 5 Newton
iterations, precision estimates (`|slope/curv|` at the converged
point) in the 1e-7 to 1e-3 substitutions/site range** - LG now
indistinguishable from the other four. This also finally answers
Lasse's original three convergence questions cleanly (previously
confounded by the LG bug):
- `MAX_ITERATIONS=50` is never approached (observed max: 5) - large
  headroom, cheap insurance, no change needed.
- `CONVERGENCE_TOL=0.1` is not too loose - achieved distance precision
  is consistently tiny (median ~1e-5, worst case ~1e-3) relative to
  typical distances (0.1-2).
- `MIN_DISTANCE`: tightened `0.01` -> `0.001` per Lasse's suggestion.
  Real data (both datasets, post-fix) never actually hits the floor at
  either value, so this isn't evidenced as a correctness fix, but
  0.001 is closer to PHYLIP protdist's own floor (0.00001) and gives
  closely-related pairs more room before clamping - low-risk
  precaution, not a bug fix.

`fastprot_mpi`'s separate `ModelMatrix.cpp` copy was NOT updated with
`mean_substitution_rate()` (unlike the earlier JTT data fix, which was
mirrored there). Decided (2026-08-11): leave it - unlike the JTT fix
(a trivial, mechanically-verifiable one-line data correction), this
one touches three files and can't be verified at all here (no MPI
available, performance work already paused there per the
2.0.0-beta.1 changelog). Deferred alongside the JTT-adjacent gaps
already tracked in "Explicitly out of scope" below.

Verified: full rebuild + `ctest` 6/6 on Clang Release and real GCC 14;
`RunExamples.sh` clean after regenerating `ex18.out` (the `-D WAG -m`
example - WAG's mean rate is 95.28, not the old flat 100, so its
distances shift by the expected small amount); `RunCliChecks.sh`
clean; direct re-check against PHYLIP's `protdist` confirmed JTT
unchanged (`2.008679`, exact match - JTT's true rate was already
~100.6, close enough that this fix barely moves it).

## Symmetrized decomposition + a build-once API (2026-08-11)

Asked whether there was more single-threaded speedup available (not
multithreading - deliberately deferred, see below), the answer was
`MatrixExpm`'s one-time decomposition itself: re-measured against the
*current* (post all prior optimizations) per-pair loop, since Q2's old
"≤0.3% of wall time, negligible" finding was measured against the much
slower pre-Eigen loop and never re-checked after the loop got 2-3x
faster. `benchmarks/bench_decomp_share.cpp`: decomposition is **92-95%
of total time at N=2 (1 pair), 55-63% at N=5, 21-26% at N=10, 6-7% at
N=25, 0.03% at N=300** - real and dominant for small datasets, still
genuinely negligible for large ones.

**Symmetrization**: every rate matrix `MatrixExpm` decomposes is
time-reversible (`eq[i]*Q(i,j) == eq[j]*Q(j,i)` by construction, true
for all 5 named models), so `Q` is similar to the symmetric matrix
`S = D^0.5*Q*D^-0.5` (`D = diag(eq)`) - decomposing `S` via
`Eigen::SelfAdjointEigenSolver` instead of the general (possibly
complex-eigenvalue) `EigenSolver` directly on `Q`, with `Q`'s own
eigenvectors/inverse recovered as `V = D^-0.5*U`, `V^-1 = U^T*D^0.5`
(no explicit matrix inverse needed - `U` is orthogonal), measured
**3.15-3.8x faster** (`benchmarks/bench_symmetric_decomp.cpp`), exact
to machine precision for the 4 exactly-reversible models
(WAG/JTT/Dayhoff/LG) and to 6.2e-5 (well inside this project's 5e-4
tolerance) for MVR - which turns out to have a small, genuine (~1.2%
relative) detailed-balance violation in its own published data, found
while verifying this change, not introduced by it. `MatrixExpm` kept
its old general-`EigenSolver` constructor too (unchanged), since
`Matrix::expm()` (the ED path, `ExpectedDistance.cpp`) calls it with
no equilibrium distribution available and wasn't in scope to touch;
the new, faster constructor takes `eq` and is what the ML path uses.

**Build-once API**: symmetrization alone wasn't the full ask - the
request was an API where decomposition is "at worst a start-up cost,"
true "regardless of whether distance estimation is done once or many
times." This wasn't hypothetical: `fastprot -b N` bootstrapping
(`fastprot/main.cpp`) calls `calculate_ml_dists()` once for the
original data plus once per replicate, and before this change *every*
call rebuilt `MatrixExpm` from scratch even though the model - and
therefore the decomposition - never changes across replicates. A
second, unbounded-call-count case was also raised directly: a future
`fastphylo-py` binding that won't know in advance how many distance
estimations will be requested, where the same guarantee still has to
hold.

Fix: `calculate_ml_dists()`'s signature changed from taking a
`model_type` (and building its own `MatrixExpm` internally, every
call) to taking an already-built `const MatrixExpm&` - it no longer
has the information needed to redo the expensive step even if a
caller wanted it to. A new `build_ml_decomposition(model_type)`
factory (`ProtDistCalc.hpp/.cpp`) does what `calculate_ml_dists()`
used to do internally, as the one explicit "start-up cost" call every
caller (CLI, tests, future bindings) makes once and reuses.
`fastprot/main.cpp`'s bootstrap loop now builds it once in `main()`,
before `processRuns()`'s dataset/bootstrap loop even starts, and
threads a `const MatrixExpm *` through to the ML branch;
`calculate_distances()`'s two dispatcher overloads (the existing
one-shot convenience path, unaffected in contract) now call
`build_ml_decomposition()` internally, immediately followed by
`calculate_ml_dists()`, preserving their exact existing "compute once"
behavior for any other caller.

**Measured end-to-end** (interleaved, `--seed` fixed for identical
resampling): on the 25-sequence globin family with `-b 300` (300 pairs/
replicate, per-pair work dominates regardless), the improvement is
small, as expected. On a 5-sequence dataset with `-b 2000` (10 pairs/
replicate - decomposition dominated before this fix): **~1.9x faster**
(0.22-0.23s -> 0.11-0.12s). On a 2-sequence dataset with `-b 5000` (1
pair/replicate, the extreme case): **~7.75x faster** (0.31s -> 0.04s).
Correctness: full rebuild + `ctest` 6/6 on Clang Release and real GCC
14; `RunExamples.sh` clean after regenerating `ex18.out` (WAG shifted
~1e-6, the same float-reordering-noise magnitude as prior `float`-path
disclosures); `RunCliChecks.sh` clean; direct byte-level comparison of
all 5 models' full 25x25 distance matrices against the pre-change
build confirmed WAG/JTT/Dayhoff/LG agree to max 1e-6 absolute (noise)
and MVR to max 6.5e-3 absolute / ~0.5% relative (the disclosed,
pre-existing detailed-balance-violation tradeoff above, propagated
through the Newton solver into the final distance - small, real, and
was flagged before implementing, not discovered after).

Deliberately not pursued in this round (asked about, explicitly
deferred): multithreading the pairwise loop - real, but a different
kind of change (thread-safety, opt-in flag design) the user wanted to
set aside in favor of single-threaded cleverness first.

## Follow-up profiling round: replacing rand() (2026-08-11)

Asked for a fresh `sample`-based profile after the decomposition fix
above, two runs: a large-N (1000 sequences) steady-state case, and a
small-N/huge-`-b` stress case built specifically to confirm the
decomposition fix under a profiler, not just wall-clock timing (it
did - no `EigenSolver` symbols appear at all in the stress-case
profile). The large-N profile found nothing further worth chasing -
~95%+ of samples are inside legitimate Eigen GEMM work or the
necessary tally-counting to build `N`. But the stress-case profile
surfaced a new, previously-hidden cost: `bootstrap_sequences()`
(`ProtSeqUtils.cpp`, pre-existing, untouched this session) calls
libc's `rand()` once per sequence-column per replicate - combined with
its own loop overhead, resampling was **~31% of total time** in that
extreme (N=2, `-b 200000`) workload. Invisible before only because the
much larger decomposition cost was masking it.

The same `seqlen * 1.0 * rand() / (RAND_MAX + 1.0)` pattern turned out
to be duplicated identically in 3 places - `ProtSeqUtils.cpp`
(fastprot, protein), `Sequence.cpp` and `Sequences2DistanceMatrix.cpp`
(both DNA, used by fastdist) - all fixed together for consistency
rather than just the one profiled.

**First attempt was a real, measured regression, not an improvement**:
`std::mt19937` + `std::uniform_int_distribution` (the obvious "modern
`<random>`" swap) measured **2.45x slower per call** than `rand()`
(17.6ns vs 7.2ns, `benchmarks/bench_bootstrap_rng.cpp`) - confirmed
end-to-end too (0.41s -> 0.57s on the `-b 50000` stress test).
`uniform_int_distribution`'s bias-avoidance machinery costs more than
`rand()`'s external-call overhead it was meant to save. Caught by
re-measuring rather than trusting the assumption that any `<random>`
engine trivially beats `rand()`.

**Root-caused and fixed properly**: the actual old code already used
its own float-scaling trick (`n * draw() / (max+1.0)`) to get a
uniform value from `rand()`'s range - the fix that actually works is
keeping that same trick, just feeding it `std::mt19937` instead of
`rand()`, skipping `uniform_int_distribution` entirely. Measured
**3.71ns/call, 1.94x faster than `rand()`**, standard-library only, no
custom code. A hand-rolled `xorshift32` measured faster still
(2.04ns/call, 3.52x) but was deliberately not adopted after being
asked for a reference - Marsaglia's base 3-shift xorshift family is
linear over GF(2) and has documented weaknesses (fails binary-rank/
linear-complexity tests), which is specifically why later variants
(xorshift+, xoshiro) add a scrambling step; standard, well-vetted
`mt19937` was preferred over a faster but weaker non-standard engine.
New shared `fastphylo/core/random_utils.hpp/.cpp`
(`seed_rng()`/`uniform_random_index()`) replaces `srand()`/`rand()` at
all 3 call sites plus their 2 seeding sites (`fastprot`/`fastdist`
`main.cpp`).

Measured end-to-end (interleaved, fixed `--seed`): a real but modest
~7-8% improvement on the extreme `-b 50000` stress case (0.41s ->
0.38s user time) - smaller than the isolated 1.94x per-call number
since the RNG draw was never the *only* cost even in that extreme
case (Eigen GEMM, `bootstrap_sequences()`'s own loop overhead, and
per-replicate sequence re-encoding all still contribute).

Verified: `ctest` 6/6 on Clang Release and real GCC 14;
`RunExamples.sh`/`RunCliChecks.sh` clean, including example 19's
bootstrap-through-`fnj` pipe (matches byte-for-byte despite the RNG
change - the fixture checks tree *topology*, not raw distances, and
this small 4-taxon dataset's phylogenetic signal is strong enough that
every bootstrap replicate converges to the same topology regardless of
which RNG produced the resampling, confirmed to be the reason, not a
coincidence); direct `--seed` reproducibility check (same `--seed`
twice gives byte-identical output, just a different resampling
sequence than the pre-change `rand()`-based code would have produced -
a real, disclosed break in cross-version reproducibility for anyone
relying on exact `--seed` output, not a bug).

## Explicitly out of scope (carried over from the investigation plan)

- `fastprot_mpi`'s separate, unsynced `Matrix`/`MaximumLikelihood`/
  `ModelMatrix` copies - performance work paused per the 2.0.0-beta.1
  changelog; also now missing PR #2's correctness fix and the LG
  mean-substitution-rate fix above (known, disclosed gaps, not
  silently forgotten - the JTT data fix *was* mirrored there earlier,
  as the one trivial, low-risk exception).
- `ExpectedDistance.cpp`'s `posterior_probability()` - a structurally
  similar allocation-churn problem, deliberately not folded in here;
  its own investigation if/when it becomes a priority.
- `count_replacements()`'s primitive - already SIMD-accelerated,
  small fraction of ML time, not revisited.
- IQ-TREE2 model-catalog sourcing, `modelestimator`-output parsing,
  `fastphylo-py` Python bindings - all real, all deliberately not part
  of this performance plan; each is its own effort when picked up.
