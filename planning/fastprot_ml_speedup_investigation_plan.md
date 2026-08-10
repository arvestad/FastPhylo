# Plan: `fastprot` maximum-likelihood speedup investigation

## Goal

Scoped for **after** 2.0.0-beta.1 ships - not blocking it, but worth
having a concrete procedure ready to start from. Direct request: "a
thorough procedure to figure this out... data-driven, with profiling
and numerical experiments," covering four specific questions (matrix
backend - MKL vs. Eigen vs. the current hand-rolled `Matrix`/LAPACK
interface; algorithmic redundancy in the spectral decomposition;
Newton-Raphson vs. a derivative-free solver like Brent's method;
input-independent precomputation) plus "there is more to consider."
This plan is a *procedure*, not a set of conclusions - every question
below gets settled by measurement against this codebase's actual data
and actual code, not by general reasoning about which library "should"
be faster.

## Current state (verified before writing this plan, not assumed)

The `speed2026a` branch (merged into `master` as part of 2.0.0-beta.1)
already did two rounds of ML-specific work - `phase0_audit.md`'s "ML
speedup round" and "part 2" - worth knowing precisely what's already
done so this investigation doesn't re-discover it:

- **Q's eigendecomposition is already cached once per
  `calculate_ml_dists()` run**, not per Newton-Raphson iteration and
  not per pair (`Matrix.hpp`'s `MatrixExpm`, constructed once in
  `ProtDistCalc.cpp::calculate_ml_dists()`, threaded through
  `likelihood_calc()`/`likelihood_deriv()`). This eliminated the
  original LAPACK-dominated profile entirely (confirmed via profiling
  at the time - no `DGEEV`/`DLAHQR`/`DTREVC` symbols in the post-fix
  profile's top 20 leaf frames). **Do not re-propose "cache the
  decomposition" as a finding - it's done.**
- `MatrixExpm::at(t)`'s diagonal-scaling step and
  `likelihood_deriv()`'s `elem_mult`+`elem_div`+`sum` chain were both
  already fused to avoid intermediate `Matrix` temporaries.
- What's left, per that round's own final note (not yet acted on):
  `likelihood_deriv()` still heap-allocates 2 `Matrix` temporaries per
  call (`pt = Qdecomp.at(t)`, `direction = Matrix::mult(pt, Q)`, both
  20x20), and every element access anywhere in `Matrix` goes through a
  bounds-checked `operator()` that throws on failure - this is
  `[[project-matrix-implementation-investigation]]`'s finding
  (2026-08-01 memory entry), not yet independently re-confirmed against
  current (post-2.0.0-beta.1) code.
- **The current "Newton-Raphson" solver is already a finite-difference
  method, not a true analytic-gradient Newton's method**:
  `likelihood_calc()` in `MaximumLikelihood.cpp` evaluates
  `likelihood_deriv()` (the score function, i.e. d/dt log-likelihood)
  at `t` and `t+delta`, finite-differences those two values to
  approximate the *second* derivative, and takes a Newton step on the
  score function's root. After the first iteration, each subsequent
  iteration needs exactly one new `likelihood_deriv()` evaluation
  (`l_new` becomes next iteration's `l_d`) - so it's already
  effectively "one function evaluation per iteration," the same
  asymptotic cost class a derivative-free method like Brent's would
  have. This matters for question 3 below: Brent's advantage here, if
  any, would have to come from *convergence robustness* (fewer
  iterations to reach tolerance, or handling the cases that currently
  trigger the method's ad hoc bail-outs - see below), not from
  "avoiding an expensive gradient," since there isn't one today.
- The solver has three hand-written bail-out conditions (`t < 1` →
  return 1, `t > 500` → return 500, `fabs(l_d) < fabs(l_new)` →
  "derivative is getting larger.. return t") beyond the normal
  tolerance check. These read as guards against divergence/
  oscillation - worth treating as evidence the current method doesn't
  always converge cleanly, and specifically checking whether a
  bracketing method would hit these same edge cases more gracefully.
- `Q` itself (the raw, un-decomposed rate matrix) is rebuilt from
  ~400 hardcoded literal doubles per model (`ModelMatrix.cpp`'s
  `get_wag()`/`get_jtt()`/etc., one function per model, called once per
  `calculate_ml_dists()` run) - plain assignment, not a computation, so
  this specific step is very unlikely to matter and doesn't need
  profiling time spent on it. What *would* be worth precomputing
  further (question 4 below) is the decomposition *result*
  (eigenvalues/eigenvectors), not the raw matrix.
- `fastprot_mpi` has its own separate, unsynchronized copies of
  `Matrix.cpp`/`MaximumLikelihood.cpp`/`ModelMatrix.cpp`
  (`src/c++/apps/fastprot_mpi/`) that never received the `speed2026a`
  fixes above (couldn't be verified - no MPI available - same standing
  limitation as the rest of this engagement). Given `fastprot_mpi` is
  now explicitly paused (2.0.0-beta.1's changelog), this investigation
  is scoped to the main library copy
  (`src/fastphylo/protein/`) only.

## The four questions, made concrete and measurable

### Q1 - Matrix backend: hand-rolled `Matrix` vs. Eigen vs. MKL

**These are two different axes, not three options on one axis - keep
them separate or the experiment design will be confused:**

- *API/allocation layer*: today's hand-rolled `Matrix` (heap-allocated
  `std::vector<double>` per instance, bounds-checked `operator()`) vs.
  Eigen's fixed-size `Eigen::Matrix<double,20,20>` (stack-allocated,
  no bounds check by default, expression templates fold chained
  operations like `elem_mult`+`elem_div`+`sum` into one pass
  automatically instead of needing it done by hand as
  `likelihood_deriv()` already does today).
- *BLAS/LAPACK vendor*: today links whatever `FindBLAS`/`FindLAPACK`
  finds (Accelerate on macOS, typically reference LAPACK or OpenBLAS
  on Linux) vs. Intel MKL as an alternative vendor for the same calls
  (`dgeev_` etc.).

Given the profiling evidence above - the *current* bottleneck is
allocation/bounds-checking overhead on tiny (20x20) matrices, not raw
FLOP throughput - the leading hypothesis is that the API/allocation
layer is the axis that matters, and the BLAS vendor is close to
irrelevant at this matrix size (vendor kernel differences show up in
FLOP-bound work; a 20x20 `dgeev_` call is dominated by fixed overhead
regardless of vendor). **Test this directly, don't assume it**: a
microbenchmark comparing (a) current `Matrix`, (b) `Matrix` with
`operator()` bounds-checks compiled out / `Eigen::Matrix<double,20,20>`
as a drop-in, (c) the same with MKL linked instead of the current BLAS
backend, run head-to-head on the actual operations
`likelihood_deriv()`/`MatrixExpm::at()` perform, at the actual 20x20
size. If the hypothesis holds, MKL adds a new, heavier build/runtime
dependency (a real cost - `INSTALL`'s dependency list, binary size,
`vcpkg`/CI provisioning on all three platforms) for a benefit this
matrix size can't realize; Eigen is header-only and already scoped as
a candidate in the pre-existing memory finding, with real prep work
listed there (confirm `Matrix`'s usage scope beyond `fastprot`,
confirm Eigen's `EigenSolver` - the general, non-symmetric case, since
these rate matrices aren't symmetric - produces results compatible
with LAPACK's `dgeev_`, byte-identical or within a documented and
justified tolerance).

### Q2 - Further reducing spectral-decomposition work

Already once-per-run, not once-per-pair or once-per-iteration (see
above) - the open question is whether "once per run" is still more
than necessary. Two sub-questions, both answerable by measurement
rather than argument:

1. **How much does that one decomposition actually cost relative to a
   typical run?** If `fastprot` is normally run on datasets large
   enough that the per-pair Newton-Raphson work dominates total time,
   amortizing one `dgeev_` call across thousands of pairs already makes
   it negligible - not worth further engineering. If `fastprot` is
   commonly run on small inputs (a handful of sequences), the one-time
   decomposition could still be a meaningful fraction of total wall
   time. Needs profiling across a realistic *range* of dataset sizes,
   not just the large benchmark datasets `speed2026a` used (small-input
   latency and large-input throughput are different questions).
2. **Could the decomposition be skipped entirely for the 6 named
   models** (WAG/JTT/DAY/ARVE/MVR/LG - each a fixed, published rate
   matrix, never data-dependent) **by embedding precomputed
   eigenvalues/eigenvectors as literal constants**, the same way `Q`
   itself is already embedded as literal doubles in `ModelMatrix.cpp`?
   This is the cleanest instance of "precompute regardless of input
   data" (question 4) - concrete, scoped to exactly 6 known matrices,
   and testable by literally computing and checking in the 6 constants,
   then comparing output byte-for-byte against the current
   runtime-`dgeev_` path before deciding whether the added source
   verbosity (6 sets of 20 eigenvalues + 400-entry eigenvector
   matrices) is worth however much time sub-question 1 finds this
   actually saves.

### Q3 - Newton-Raphson vs. a derivative-free solver (Brent's method)

Given the current solver is already a one-evaluation-per-iteration
finite-difference method (see above), the fair experiment is not
"gradient-based vs. gradient-free" (there's no analytic gradient being
computed today either way) but a head-to-head on the same root-finding
problem (`likelihood_deriv(t) = 0`):

- Implement Brent's method (or Brent-Dekker; needs an initial bracket
  containing the root - `kimura_distance()`'s existing starting value
  plus the solver's own `[1, 500]` clamp range are natural bracket
  candidates, worth checking they actually bracket a root in practice
  before assuming it) against the same `likelihood_deriv()` function.
- Compare, across a real sample of `N` matrices (both the small fixture
  data this project already has and, ideally, larger/messier real
  datasets): number of `likelihood_deriv()` evaluations to reach the
  same tolerance, wall time, and - specifically - behavior on whatever
  inputs currently trigger the three ad hoc bail-out conditions (does
  Brent's method converge cleanly there instead of hitting a
  heuristic early return?).
- A solver swap changes numerical results in general (different
  iteration path, possibly different converged value under finite
  tolerance) - this is **not** expected to be byte-identical the way
  every pure-performance change in this engagement has been. Define
  and document an explicit acceptable-difference tolerance *before*
  running the comparison, and report actual observed differences
  against it - don't discover the tolerance question after the fact.

### Q4 - Other input-independent precomputation

Q2's embedded-eigendecomposition idea is the concrete instance found
so far. Systematically check for others while doing this work:
anything in the per-run or per-pair path that depends only on which
`model_type` was selected (not on the actual sequences) is a
precomputation candidate. Worth a deliberate sweep of
`ProtDistCalc.cpp`'s `calculate_ml_dists()`/`calculate_ed_dists()` and
`ModelMatrix.cpp` for anything else shaped like this - this plan
doesn't enumerate more candidates now because they should come from
reading the current code freshly when this work actually starts, not
from guessing ahead of time.

## Methodology

- **Profile first, on current code, before forming hypotheses** - the
  single sharpest lesson from `phase0_audit.md`'s own record (a profile
  taken a few commits earlier gave an actively misleading picture of
  where time was going later). Re-profile `fastprot`'s ML path fresh at
  the start of this work, even though the "current state" section above
  already describes what the last profiling round found - confirm it
  still holds on 2.0.0-beta.1's actual code before building on it.
- Same profiling approach as `phase0_audit.md` used (check that doc for
  the exact tool/flags - sampling profiler, optimized build, realistic
  dataset) for continuity/comparability with the existing record.
- **Benchmark across a range of dataset sizes**, not just one - Q2
  specifically needs this (small-input latency vs. large-input
  throughput are different questions), and it's good practice for Q1/Q3
  too. `benchmarks/gen_dataset.py`/`run_benchmarks.sh` already exist
  from `speed2026a`'s work - extend rather than replace.
- **Verify correctness at whatever tolerance is appropriate per
  change**: byte-identical for pure refactors/allocation changes (Q1's
  Matrix-layer swap, Q2's decomposition-caching extensions), an
  explicitly documented and justified numerical tolerance for anything
  that changes the actual computation path (Q3's solver swap, Q2's
  embedded-constants approach if the embedded values are computed by a
  different code path than the current runtime one).
- Interleave old/new binary runs when measuring wall time (established
  practice from `phase0_audit.md`, given thermal-throttling drift was
  observed there) rather than trusting a single before/after pair.

## Phases (draft - refine once this work actually starts)

1. **Fresh profiling baseline** on current (2.0.0-beta.1) code, across
   a range of dataset sizes. Confirms or corrects the "current state"
   section above; produces the actual numbers Q1-Q4 get measured
   against.
2. **Q1 - matrix backend microbenchmark**: current `Matrix` vs. Eigen
   vs. (if the microbenchmark's own results justify spending the setup
   effort) MKL, on the actual 20x20 operations involved. Decision point:
   is the allocation/bounds-check-layer hypothesis confirmed? If yes,
   Eigen is the likely direction (cheaper dependency than MKL, matches
   the pre-existing memory finding); if the BLAS vendor turns out to
   matter more than expected, revisit.
3. **Q2 - decomposition-cost-vs-dataset-size sweep**, then (if
   justified by that data) the embedded-constants experiment for the 6
   named models.
4. **Q3 - Brent's-method experiment**, run independently of Phases 2/3
   (touches the solver, not the matrix layer - low interference risk,
   could happen in parallel).
5. **Synthesis**: given Phases 1-4's actual measurements, write the
   real implementation plan (which changes are worth making, in what
   order, what each one's expected impact is) - this document is the
   investigation procedure, not that plan.

## Explicitly out of scope

- `fastprot_mpi`'s separate `Matrix`/`MaximumLikelihood`/`ModelMatrix`
  copies - paused, per 2.0.0-beta.1's changelog; not touched by this
  investigation regardless of what it finds for the main library.
- `ExpectedDistance.cpp`'s `posterior_probability()` - flagged in the
  same memory entry as having a "structurally identical allocation-
  churn problem," explicitly deprioritized in the `speed2026a` round
  that found it. Revisit only if Q1's Matrix-layer work ends up
  touching shared code `ExpectedDistance.cpp` also uses, in which case
  it may improve "for free" - not a goal to pursue independently here.
- `count_replacements()`'s own primitive - already SIMD-accelerated and
  profiled as a small fraction of total ML time in the current
  pipeline; not expected to be where further gains live, per the
  "current state" section above.
