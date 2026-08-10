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

**Updated 2026-08-10** with direct follow-up feedback: which dataset
size is actually the hard case (reprioritizes the four questions - see
"Scope decisions" below), a `float`-vs-`double` question folded into
Q1, a concrete numerical tolerance for Q1/Q3's non-byte-identical
comparisons, and four related-but-distinct scope items (removing the
`ARVE` model, sourcing more models from IQ-TREE2's catalog, reading
`modelestimator` output as a custom model, and making `fastphylo-py`
benefit from this work) that inform Q2/Q4 without being performance
questions themselves.

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

## Scope decisions from direct feedback (2026-08-10), before the four questions

- **The hard case is large `N` - hundreds to thousands of sequences,
  not a handful.** This directly reprioritizes the four questions:
  anything that costs *per pair* or *per Newton-Raphson iteration*
  (Q1's allocation/bounds-check overhead, Q3's solver choice) scales
  with `N²`-ish and dominates at this size; anything that's a fixed
  *per-run* cost (Q2's one-time spectral decomposition) gets amortized
  across the ~`N²/2` pairs in a large run almost regardless of how
  expensive it is. **Q1 and Q3 are now the primary levers for the
  stated hard case; Q2 is very likely a secondary win** (real, and
  still worth doing per the feedback below, but shouldn't be expected
  to move the needle on large-`N` runs the way Q1/Q3 could). Profiling
  (Phase 1) should prioritize realistic large-`N` datasets accordingly,
  not just extend the size range for its own sake.
- **Remove the `ARVE` model** ("non-standard"). **Done** (this branch):
  `get_arvestad()`/`get_arve_eq()` and the `arve` enum case deleted
  from both the main library (`ModelMatrix.hpp`/`.cpp`,
  `ProtDistCalc.cpp`) and `fastprot_mpi`'s separate, otherwise-
  unsynced copy (kept in step for this one mechanical change, unlike
  the performance work, since it's low-risk and just textual);
  `fastprot`'s and `fastprot_mpi`'s `-D`/`--distance-function` CLI help
  text and `CheckedTransformer` maps updated to match. Verified: full
  rebuild (Clang + real GCC 14) + `ctest` 4/4 + `RunExamples.sh` 20/20
  byte-identical (no fixture used `ARVE`) + `RunCliChecks.sh`, a
  `-DWITH_LIBXML=OFF` build, and a direct check that `-D ARVE` is now
  rejected by CLI11's `CheckedTransformer`. Shrinks the named-model set
  Q2/Q4's precomputation ideas apply to from 6 to 5 (`JTT`/`DAY`/`MVR`/
  `WAG`/`LG`).
- **Add more matrices, using IQ-TREE2's model catalog as the
  reference.** A real, separate sourcing task - IQ-TREE2 supports
  many more protein models (e.g. Dayhoff-DCMut, VT, PMB, Blosum62, the
  mtREV/mtMAM/mtART family, HIVb/HIVw, FLU, rtREV, cpREV, the
  Q.pfam/Q.bird/Q.insect/... family, LG4M/LG4X, and more) - and this
  plan deliberately does not attempt to enumerate or transcribe exact
  published rate values here from memory; that needs to be sourced
  carefully (IQ-TREE2's own documentation/source, or the original
  papers) and verified the same way this engagement has verified
  everything else, not guessed. Directly relevant to Q2/Q4: more named
  models means a bigger, more valuable precomputed-decomposition table
  if that direction pans out - but the sourcing work itself is outside
  this document's scope and deserves its own plan when picked up.
- **Support reading `modelestimator` output as a custom (user-supplied)
  model.** A real feature request, and it sharpens Q2/Q4's scope: a
  user-supplied rate matrix is inherently data-dependent (read from a
  file at runtime, potentially different on every invocation) and can
  **never** be a precomputation candidate the way the fixed named
  models are - it still benefits from every other speedup this
  investigation finds (Q1's allocation work, Q3's solver, the
  already-done once-per-run decomposition caching), just not from
  Q2/Q4's "skip the decomposition entirely" idea. Worth keeping this
  distinction explicit wherever Q2/Q4 gets implemented, so "precompute
  for known models" code doesn't accidentally get asked to handle an
  arbitrary custom matrix too.
  **Found while doing the `ARVE` removal above**: `fastprot`/
  `fastprot_mpi` already have a `-F`/`--model-file` option
  ("Read matrix and equilibrium distribution from file, when used
  --distance-function is disregarded") - a custom-model mechanism may
  already exist in some form. Not yet checked whether its file format
  matches (or could be made to match) `modelestimator`'s actual output
  - that's the concrete next step for this item, not building the
  custom-model path from scratch.
- **`fastphylo-py` should benefit from this work.** Investigated
  2026-08-10 (public at `github.com/arvestad/fastphylo-py`, pybind11 +
  scikit-build-core, ships a compiled `_fastphylo` extension). Findings
  that change this plan's shape:
  - **DNA/tree side**: `src/cpp/fastphylo/` is a **vendored, stale copy**
    of this library's core - it still has `Object.hpp`/`.cpp` and
    `Exception.hpp`/`.cpp`, both of which this repo deleted (the
    Object-modernization and legacy-error-handling work). It is not a
    submodule, dependency, or otherwise synced to this repo - it was
    forked at some pre-modernization point and has diverged since. It
    is **not** a live consumer of anything this plan speeds up unless
    someone re-syncs it or switches it to link the real library.
  - **Protein ML side is a from-scratch reimplementation, not a wrapper
    around `fastprot`'s code at all** - there is no `Matrix.hpp`,
    `ModelMatrix.hpp`, or `MaximumLikelihood.cpp` in `fastphylo-py`.
    Instead, `src/fastphylo/protein.py` defines the rate matrices and
    does the **eigendecomposition in Python** (numpy/scipy), caching
    each named model's decomposition once. Only the branch-length
    optimization is in C++: `compute_protein_distances_cpp(names, seqs,
    right_e, left_e, evals)`, taking the already-decomposed
    eigenvectors/eigenvalues as plain numpy arrays and running
    **Brent's bounded minimizer** per pair.
  - **This bears directly on Q3**: `fastphylo-py` has already, on its
    own, replaced this repo's finite-difference "Newton-Raphson" with
    Brent's method for the same branch-length-optimization problem.
    That's independent evidence (not just the user's hunch) that Q3's
    experiment is worth doing, and a working reference implementation
    to compare correctness/speed against instead of building the
    Brent's-method experiment from a blank page.
  - **This bears directly on Q1 and the "more models"/`modelestimator`
    item**: `fastphylo-py`'s Python side already carries 16 named
    protein models (WAG, JTT, JTT-DCMUT, LG, VT, HIVb, HIVw, cpREV,
    BLOSUM62, Dayhoff, DCMUT, MtREV, RtREV, PMB, + 2 more) - more than
    this repo's 9 post-`ARVE`-removal - and already has the exact split
    this plan was trying to invent for `modelestimator` support:
    *decomposition* (from a named model, or in principle from any
    externally-supplied matrix) is separate from *optimization*, and
    only the decomposition's output (eigenvectors + eigenvalues)
    crosses the Python/C++ boundary. That argues for shaping `fastprot`
    itself the same way - a thin, fast optimizer core that consumes a
    precomputed decomposition, fed either by the library's own built-in
    models or by an external one (`-F`/`--model-file`,
    `modelestimator`) - rather than keeping decomposition and
    optimization welded together as they are in today's
    `MaximumLikelihood.cpp`. If that reshaping happens, the *same* C++
    optimizer core becomes a natural candidate for `fastphylo-py` to
    link against directly (replacing its standalone
    `compute_protein_distances_cpp`), which is the actual "benefit from
    this work" the user asked for - not exposing `Matrix`/Eigen to
    Python at all, just the decomposition-in/optimize-out boundary
    `fastphylo-py` already assumes.
  - Not yet investigated, deliberately out of scope for Phase 0: whether
    `fastphylo-py`'s Python-side eigendecomposition (numpy/scipy) is
    itself a bottleneck at the "hundreds to thousands of sequences"
    scale, and whether re-syncing or replacing its vendored C++ core is
    worth doing at all - both are calls to make once Q1/Q3 have real
    numbers, not now.

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

**Third axis, direct question raised in feedback: does `float` instead
of `double` help?** Plausible, worth testing, folded into this same
experiment rather than treated separately - it's really "what scalar
type does the matrix backend use," the same shape of question as (a)
vs. (b) above. The case for it:

- The tolerance guidance from Q3 below (agreement to 3 decimal places
  is already more precision than the data can justify) applies just as
  well here - `double`'s ~15-16 significant digits is far beyond what
  3 decimal places of a distance value need, so `float`'s ~7
  significant digits has a lot of headroom before it could plausibly
  matter, *if* error doesn't get amplified somewhere in the pipeline.
- `float` halves memory per `Matrix` (20x20: 1600 bytes vs. 3200) and
  doubles the number of lanes available per SIMD register, which could
  matter for `likelihood_deriv()`'s per-element arithmetic loop if it
  vectorizes (auto-vectorized or via Eigen, which vectorizes chained
  expressions by default) - though at 20x20 both sizes fit comfortably
  in L1 cache, so cache-footprint effects specifically are likely
  negligible; the SIMD-width effect is the more plausible win.

The case for caution, and why this needs measurement rather than
assumption:

- LAPACK's `dgeev_` is double-only; there's a separate single-precision
  routine (`sgeev_`), meaning a `float` path either needs that separate
  routine (a real, if small, extra branch to verify) or needs the
  decomposition to go through Eigen's own `EigenSolver` instead of
  LAPACK at all - tying this axis to the Eigen decision in (a)/(b)
  above rather than making it fully independent.
- Eigendecomposition of a matrix that then feeds a matrix *exponential*
  is exactly the kind of numerically-sensitive operation where small
  input errors don't necessarily propagate linearly - "3 decimal places
  is enough for the final answer" doesn't automatically mean "any
  precision loss mid-pipeline is fine," and needs checking empirically
  against real `Q` matrices and real `t` ranges, not assumed from the
  end-to-end tolerance alone.

Concretely: once Eigen is in the microbenchmark for (a)/(b), templating
the same benchmark on `float` as well as `double` is close to free
(Eigen matrices are templated on scalar type already) - do both, and
report accuracy-vs-speed for `float` the same way (a)/(b)/(c) get
reported, using the 3-decimal-place tolerance from Q3 as the
correctness bar.

### Q2 - Further reducing spectral-decomposition work

Already once-per-run, not once-per-pair or once-per-iteration (see
above) - the open question is whether "once per run" is still more
than necessary. **Given the hard case is large `N` (see "Scope
decisions" above), expect this question's answer to be "already
negligible" more often than not** - one `dgeev_` call amortized across
the ~`N²/2` pairs a large run does is unlikely to be where the time
goes. Still worth answering with real numbers rather than skipping
(explicitly requested: "an interesting opportunity for speedup"), and
it directly helps the *small*-`N` case even if large-`N` turns out
unaffected. Two sub-questions, both answerable by measurement rather
than argument:

1. **How much does that one decomposition actually cost relative to a
   typical run, across the realistic size range** - specifically
   including the small end, since that's where this sub-question's
   answer is most likely to actually matter, per the reprioritization
   above.
2. **Could the decomposition be skipped entirely for the named models**
   (5 after `ARVE`'s removal - `JTT`/`DAY`/`MVR`/`WAG`/`LG`, more if
   the IQ-TREE2 catalog-expansion work lands first - each a fixed,
   published rate matrix, never data-dependent) **by embedding
   precomputed eigenvalues/eigenvectors as literal constants**, the
   same way `Q` itself is already embedded as literal doubles in
   `ModelMatrix.cpp`? This is the cleanest instance of "precompute
   regardless of input data" (question 4) - concrete, scoped to
   exactly the known named matrices, and testable by literally
   computing and checking in the constants, then comparing output
   byte-for-byte against the current runtime-`dgeev_` path before
   deciding whether the added source verbosity is worth however much
   time sub-question 1 finds this actually saves. **Only applies to the
   fixed named models** - a `modelestimator`-supplied custom matrix
   (see "Scope decisions" above) is inherently data-dependent and must
   always be decomposed at runtime; keep that distinction explicit in
   whatever code implements this.

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
  **Working tolerance, per direct guidance**: agreement to 3 decimal
  places is already more precision than the underlying data justifies
  - "even three decimals is a bit fishy" given real biological data's
  own noise floor. Treat two computed distances as equivalent if they
  round to the same value at 3 decimal places (pin down the exact
  comparison - absolute vs. relative difference, how ties/rounding
  edge cases are handled - when the actual experiment is set up, this
  is the guiding principle, not a fully-specified epsilon yet). This is
  also the working tolerance for Q1's `float`-vs-`double` experiment
  above - same underlying justification, same bar.

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
  `Matrix`-layer swap *at `double` precision*, Q2's decomposition-
  caching extensions), the 3-decimal-place working tolerance from Q3
  for anything that changes the actual computation path (Q3's solver
  swap, Q1's `float`-precision variant, Q2's embedded-constants
  approach if the embedded values are computed by a different code
  path than the current runtime one).
- Interleave old/new binary runs when measuring wall time (established
  practice from `phase0_audit.md`, given thermal-throttling drift was
  observed there) rather than trusting a single before/after pair.

## Phase 1 results (2026-08-10): fresh ML profiling baseline

Environment: Apple Silicon (arm64), AppleClang 15, dedicated
`build-ml-profile/` (RelWithDebInfo, `-O2 -g`, not the main `build/`
dir) so this round doesn't disturb the existing build. Synthetic
datasets via `benchmarks/gen_dataset.py`, weighted toward the stated
hard case: 100/300/600/1000 protein sequences x 300 residues
(`benchmarks/data/ml_bench_*x300.fasta`), `fastprot -D WAG -m`.

**Wall-clock scaling** (median of 3 reps, `/dev/null` output to
isolate compute from I/O):

| n_seqs | pairs | median wall time |
|---|---|---|
| 100 | 4,950 | 0.150s |
| 300 | 44,850 | 1.046s |
| 600 | 179,700 | 4.099s |
| 1000 | 499,500 | 11.092s |

Scales essentially linearly with pair count (e.g. 600→1000 is 2.78x
the pairs for 2.71x the time) - no surprises here, confirms per-pair
cost is roughly constant as `N` grows, as expected once decomposition
is cached once per run rather than per pair.

**Profiling** (macOS `sample`, 1ms interval, 1000x300 dataset, 15s
capture): confirms the once-per-run decomposition caching from
`speed2026a`'s "ML speedup round" is still doing its job - zero
`DGEEV`/`DSYEV`/`DGETRF`/`DGETRI`/`DLAHQR`/`DTREVC` samples anywhere in
the profile. But it surfaces a **new, more specific finding** than
`phase0_audit.md`'s post-caching round did:

- Flat (top-of-stack) attribution: `fastprot`'s own code 36.1%,
  `libBLAS.dylib` (Accelerate) 45.1%, `libsystem_malloc` 8.8%.
- The call graph traces that 45.1% to exactly one call site:
  `likelihood_calc` → `likelihood_deriv` → `MatrixExpm::at(t)` →
  `Matrix::mult(left, eigenvectors_inv)` (`Matrix.cpp:318`, the
  "two matrix multiplies" comment left after `speed2026a`'s "part 2"
  round replaced the diagonal-factor multiply with a column-scale loop
  - this is the *other*, still-dense, multiply) → `dgemm_` → deep
  inside Accelerate's `DGEMM`/dispatch machinery. Within
  `likelihood_deriv`'s own samples, 65.6% are inside `MatrixExpm::at`;
  within that, 92.8% are inside `Matrix::mult`; within *that*, 99.5%
  are inside `DGEMM` itself, not call/dispatch overhead layered
  outside it.
- **This is a 20x20 matrix multiply** - 8,000 multiply-adds, trivial
  for even naive scalar code - going through a general-purpose BLAS
  implementation via Accelerate that's built for problems orders of
  magnitude larger. It's the single largest attributable cost in the
  whole ML pipeline today, ahead of `likelihood_deriv`'s own
  allocation/`elem_mult`/`elem_div` arithmetic that `speed2026a`'s
  "part 2" round already reduced from 6 `Matrix` temporaries to 3.

**What this changes for Q1**: the plan's existing Q1 framing (hand-
rolled `Matrix` vs. Eigen, allocation/bounds-check overhead as the
hypothesized culprit) is still right in direction but was aimed at the
wrong primary target - the allocation overhead is real (`malloc` is
8.8%) but smaller than this one `dgemm_` call site. A fixed-size
`Eigen::Matrix<double,20,20>` multiply (compile-time unrolled, no
external call, no BLAS dispatch/thread-pool/blocking logic sized for
large matrices) is a strong candidate to eliminate most of that 45.1%
directly, independent of whatever the allocation-layer swap saves
separately. Phase 2's microbenchmark should measure this specific
operation (20x20 dense multiply, called at the same rate `MatrixExpm::
at` calls it) explicitly, not just generic elementwise ops.

Not yet done, left for Phase 2: confirming Eigen (or a hand-rolled
loop) is actually faster here in practice rather than just in theory -
tiny-matrix BLAS overhead is a known general pattern, but this
codebase's actual numbers are what settles it, per this plan's own
"measurement, not general reasoning" rule.

## Phases (draft - refine once this work actually starts)

0. **Scope-settling prerequisites**, independent of profiling: remove
   `ARVE` (small, self-contained) - **done**; investigate `fastphylo-py`
   (where it lives, current architecture, what consuming `libfastphylo`
   would take) - **done**, see "Scope decisions from direct feedback"
   above. The IQ-TREE2 model-catalog expansion is *not* included here -
   real sourcing work, its own plan when picked up, not a prerequisite
   for the performance questions. Phase 0 is complete; ready to start
   Phase 1.
1. **Fresh profiling baseline** on current (2.0.0-beta.1) code, across
   a range of dataset sizes **weighted toward the stated hard case**
   (hundreds to thousands of sequences), not evenly spread. Confirms or
   corrects the "current state" section above; produces the actual
   numbers Q1-Q4 get measured against. **Done** - see "Phase 1 results"
   above; found the per-iteration 20x20 `dgemm_` call in `MatrixExpm::
   at` is now the single largest cost (45.1% of samples), sharpening
   Q1's target beyond the allocation-layer hypothesis.
2. **Q1 - matrix backend microbenchmark** (highest expected impact,
   given the hard case is large `N`): current `Matrix` vs. Eigen vs.
   (if the microbenchmark's own results justify spending the setup
   effort) MKL, each at both `double` and `float`, on the actual 20x20
   operations involved. Decision point: is the allocation/bounds-check-
   layer hypothesis confirmed? If yes, Eigen is the likely direction
   (cheaper dependency than MKL, matches the pre-existing memory
   finding); does `float` clear the 3-decimal-place tolerance bar, and
   if so does it add a meaningful further win on top of the backend
   choice, or is the allocation-layer fix where all the gain already
   was? If the BLAS vendor turns out to matter more than expected,
   revisit.
3. **Q3 - Brent's-method experiment** (also high expected impact, same
   reason as Q1 - scales with pairs × iterations). Independent of
   Phase 2 (touches the solver, not the matrix layer) - could run in
   parallel with it. `fastphylo-py`'s `compute_protein_distances_cpp`
   (see "Scope decisions" above) already does exactly this swap for the
   same problem - use it as a working reference for correctness/speed
   comparison rather than designing the experiment from scratch.
4. **Q2 - decomposition-cost-vs-dataset-size sweep**, then (if
   justified by that data) the embedded-constants experiment for the
   named models. Lower expected impact at the stated hard case than
   Phases 2-3 (see "Scope decisions" above) but cheap to answer and
   worth having the real numbers for, plus it directly helps the small-
   `N` case.
5. **Synthesis**: given Phases 1-4's actual measurements, write the
   real implementation plan (which changes are worth making, in what
   order, what each one's expected impact is) - this document is the
   investigation procedure, not that plan. Also the point to revisit
   `fastphylo-py`'s findings from Phase 0 in light of whatever Q1
   concluded, if bindings turn out to be relevant to how any of this
   actually ships.

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
- Actually sourcing the IQ-TREE2 model catalog, and actually
  implementing `modelestimator`-output parsing - both real, scoped as
  their own future work in "Scope decisions" above, not detailed here.
- Building `fastphylo-py` Python bindings, or any other concrete
  integration work - Phase 0 only investigates and reports; committing
  to a specific integration approach is a decision for after that
  investigation, not part of this plan.
