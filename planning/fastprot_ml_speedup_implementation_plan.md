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

## Explicitly out of scope (carried over from the investigation plan)

- `fastprot_mpi`'s separate, unsynced `Matrix`/`MaximumLikelihood`/
  `ModelMatrix` copies - performance work paused per the 2.0.0-beta.1
  changelog; also now missing PR #2's correctness fix (a known,
  disclosed gap, not silently forgotten).
- `ExpectedDistance.cpp`'s `posterior_probability()` - a structurally
  similar allocation-churn problem, deliberately not folded in here;
  its own investigation if/when it becomes a priority.
- `count_replacements()`'s primitive - already SIMD-accelerated,
  small fraction of ML time, not revisited.
- IQ-TREE2 model-catalog sourcing, `modelestimator`-output parsing,
  `fastphylo-py` Python bindings - all real, all deliberately not part
  of this performance plan; each is its own effort when picked up.
