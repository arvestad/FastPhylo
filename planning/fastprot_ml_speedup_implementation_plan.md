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
   before committing to it).
2. `double` Eigen swap (steps 2-5 above) - the actual work, verified
   change, real measured win.
3. Re-profile, confirm.
4. Decide on `float` with real numbers in hand.
5. Update dependency docs, changelog entry (`2.0.0-beta.4` or
   whatever's next, per `RELEASING.md`'s established prerelease
   convention).

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
