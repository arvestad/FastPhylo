# Plan: Speed up FastPhylo protein distance computation

## Goal

Replace the internal representation of protein sequences in FastPhylo /
fastphylo-py with a SIMD-friendly, byte-per-residue encoding, use it to
accelerate (a) pairwise mismatch/Hamming-style comparison and (b) tallying
of a 20x20(+ambiguity) amino-acid replacement/substitution count matrix,
and prove — with a rigorous benchmark and a correctness harness — that the
new code is both correct and faster than the current implementation.

This document is the brief for Claude Code. It assumes Claude Code has
read/write access to the FastPhylo and fastphylo-py repositories and will
run all commands itself. Concrete file paths and function names are
intentionally left as placeholders — Phase 0 is dedicated to filling
them in from the real codebase before any code changes are made.

## Target and scope (updated after initial audit)

A first audit pass found three separate implementations of protein
pairwise comparison across the two repos, not one:

1. **FastPhylo's `fastprot` CLI tool** (`src/c++/programs/fastprot/`:
   `ProtSeqUtils.cpp`, `ProtDistCalc.cpp`, `ExpectedDistance.cpp`, plus
   `fastprot_mpi/` copies). Sequences are `Sequence` objects (`std::string
   seq`, `Sequence.hpp`). Comparison primitives: `count_id_dist()`
   (Hamming/identity fraction, feeds the id/JC/Kimura/Storm-Sonnhammer
   closed-form distance corrections) and `count_replacements()` (tallies
   a 20x20 `Matrix`, feeds the expected-distance and ML models).
2. **fastphylo-py's pybind11 bindings** (`fastphylo-py/src/cpp/
   bindings.cpp`): an independently written reimplementation
   (`protein_count_matrix()` + its own Brent-minimizer ML solver) that
   does *not* call into `fastprot`'s code, despite vendoring copies of
   some other core files (`Sequence`, `NeighborJoining`, `DNA_b128`, ...).
3. **fastphylo-py's pure-Python path** (`fastphylo-py/src/fastphylo/
   protein.py`): a third reimplementation of the count matrix and Kimura
   distance that loops per-residue in Python/NumPy.

**Decision: this plan targets #1 only — the `fastprot` CLI tool in the
FastPhylo repo, primarily its Kimura distance estimator (closed-form,
driven by `count_id_dist()`), and secondarily `count_replacements()`
since the ED/ML models share it.**

fastphylo-py (#2 and #3) is explicitly **out of scope** for this round:
- Lasse's actual ask is a faster `fastprot`, especially the Kimura path.
- fastphylo-py was believed to reuse FastPhylo's code directly; the audit
  shows it doesn't (see finding above). Making that true — turning
  FastPhylo into a proper library that fastphylo-py links against instead
  of vendoring/reimplementing — is a valuable follow-up project, but a
  separate one. Do not attempt it as part of this plan.
- Do not modify fastphylo-py's `bindings.cpp` or `protein.py` here. If
  Phase 0 profiling shows the Kimura/`count_id_dist()` path in `fastprot`
  is already a small fraction of runtime, stop and report that rather
  than expanding scope into fastphylo-py to chase a bigger number.

Everywhere below that says "FastPhylo core... exposed through the Python
bindings in fastphylo-py," read it as scoped to `fastprot`/`fastprot_mpi`
and the core types they use (`Sequence`, `Matrix`, `DistanceMatrix`,
etc.) within the FastPhylo repo only. The Phase 3 bullet about updating
fastphylo-py bindings is dropped (see Phase 3 below).

**SIMD dispatch: nothing existing to reuse.** FastPhylo's only current
SIMD code (`DNA_b128`, used by `fastdist`) is not runtime-dispatched —
it's a hardcoded `-msse2` compile flag (`src/c++/CMakeLists.txt`), with
simde translating those calls to NEON on ARM at *compile* time (see
commit `b5cef42`). There is no runtime CPU-feature detection anywhere in
FastPhylo to reuse for the new protein code; if the design calls for
runtime dispatch, that's new infrastructure to build, not a
search-and-reuse task. Separately, `fastdist`'s DNA SIMD path is old
(SSE2 only, no AVX2/AVX-512) and its real-world speedup is unclear.
Modernizing it is a plausible future project Lasse is aware of, but it
is explicitly not part of this plan — don't fold it in opportunistically.

**Ambiguity codes are currently just dropped.** `getAAInd()`
(`ProtSeqUtils.cpp`) maps only the 20 canonical amino acid letters (upper
or lower, via `toupper()`); anything else (gaps, X/B/Z, whitespace, ...)
returns the sentinel `100` and is silently excluded from both
`count_id_dist()`'s numerator/denominator and `count_replacements()`'s
tally. Phase 1 should expect to preserve exactly this "unknown symbol
contributes nothing" behavior, not invent new semantics for ambiguity
codes — the current code doesn't have any worth preserving beyond that.

**No CI and no automated test runner in this repo today.** FastPhylo has
no `.github/workflows` and `code_tests/*.cpp` isn't wired into
`CMakeLists.txt` (there's no target that builds/runs them). The closest
thing to an existing regression suite is `examples/RunExamples.sh`, which
runs the CLI binaries (`fastdist`, `fastprot`, `fnj`) and diffs their
output against `examples/expected_output/`. Phase 0/4 should account for
this: "wire into the existing test suite/CI" means extending
`RunExamples.sh`-style checks and/or getting `code_tests/` actually built
and run, not assuming a CI pipeline exists to hook into.

## Design decisions to carry into the implementation

These were settled in discussion before this plan was written; Claude
Code should not re-litigate them without a strong reason, but should
flag it clearly if the real code makes one of them infeasible.

1. **Encoding: one byte per residue (`uint8_t`), not sub-byte bit-packing.**
   Map each amino acid (plus gap/indel and standard ambiguity codes, e.g.
   X/B/Z) to a small integer code, 0..N-1 with N <= ~25. Do **not** pack
   5-bit codes into 64-bit words — symbol boundaries then don't align to
   byte boundaries, which prevents using SIMD byte-compare instructions
   directly and forces an unpack step that erases the benefit. A flat
   `uint8_t` array per sequence is simple, cache-friendly (sequences are
   typically small enough to sit in L1/L2 already), and lets SIMD compare
   instructions operate on 16/32/64 residues per instruction.
2. **Comparison primitive:** byte-wise SIMD compare (e.g. `pcmpeqb`
   equivalents via compiler intrinsics or auto-vectorization-friendly C)
   producing a bitmask, then popcount, to get mismatch/match counts across
   a full sequence pair in a handful of instructions per 32/64 residues.
3. **Replacement/substitution count matrix:** this is a scatter-increment
   into a small (≤25x25) table and does not vectorize as cleanly as the
   compare case. Default to a well-unrolled scalar loop over `uint8_t`
   codes (`counts[a * N + b]++`), since the table fits in L1 and this is
   usually not the bottleneck; only reach for per-lane histogram tricks or
   AVX-512 scatter if profiling shows this loop actually dominates runtime.
4. **No behavior change to public results.** The new code must produce
   numerically identical (or, where floating point is involved,
   equivalent within a justified tolerance) distances to the current
   implementation. This is a performance refactor, not an algorithm
   change.
5. **Old and new implementations coexist until validated.** Do not delete
   or rewrite the existing protein code path in place. Add the new code
   alongside it, behind a compile-time or runtime switch, so the
   correctness/benchmark harness can run both against the same inputs and
   compare. Remove the old path only after Phase 5 passes and you (Lasse)
   have signed off.

## Phase 0 — Audit the existing code

Do this before writing any new code. Confirm/refine the findings already
captured under "Target and scope" above rather than starting from
scratch — but verify them against the actual code before relying on
them, since that section is a summary, not gospel.

- Confirm the protein sequence representation used by `fastprot`:
  `Sequence` (`src/c++/Sequence.hpp`, `std::string seq`), shared with DNA
  (same base type, no protein-specific subclass) — validation/alphabet
  handling happens in `getAAInd()` (`ProtSeqUtils.cpp`), not at load time.
- Confirm the comparison primitives and their call sites:
  `count_id_dist()` (Hamming/identity fraction; called from
  `calculate_id_dists`, `calculate_jc_dists`, `calculate_kimura_dists`,
  `calculate_stormsonnhammer_dists` in `ProtDistCalc.cpp`) and
  `count_replacements()` (20x20 `Matrix` tally; called from
  `calculate_ed_dists`, `calculate_ed_dists_with_sd`, `calculate_ml_dists`
  in `ProtDistCalc.cpp`, feeding `ExpectedDistance.cpp`/ML code).
- **Profile a realistic `fastprot` run** (e.g. `perf`/`gprof`/`valgrind
  --tool=callgrind` on a few hundred sequences of typical length, default
  model = Kimura) to get a flat/call-graph profile before committing
  further phases of work. In particular, confirm what fraction of total
  wall time is actually spent in `count_id_dist()` for the Kimura path
  (it's a single scalar pass with a closed-form correction afterward, so
  the primitive itself may dominate — but confirm rather than assume),
  versus how much the ED/ML models spend in `count_replacements()` vs.
  the numerical routines that consume its output (`expected_distance()`,
  `likelihood_calc()`). Report this honestly even if it's a mixed answer
  (e.g. "Kimura is dominated by the primitive, ML is not").
- Identify existing tests/benchmarks for `fastprot` specifically:
  `examples/RunExamples.sh` (examples 6/7, diffed against
  `examples/expected_output/ex6.out`/`ex7.out`) and any relevant files
  under `src/c++/code_tests/` (currently unbuilt — see scope note above).
  New tests should follow these conventions, extended as needed, rather
  than inventing a parallel style or reaching for fastphylo-py's pytest
  suite (out of scope here).
- Confirm the build system for `fastprot`/`fastprot_mpi`
  (`src/c++/CMakeLists.txt`) and how to build and run it locally,
  independent of fastphylo-py.
- Confirm compilation environment for `fastprot`: current flags
  (`ARCH_FLAGS = -msse2`, see `src/c++/CMakeLists.txt`), and that (per
  the scope note above) there is no existing runtime SIMD dispatch to
  reuse — any dispatch mechanism for the new protein code is new work.

Deliverable for this phase: a short written summary (a few paragraphs,
committed as a markdown note or PR description) of what was found —
especially the profiling results — so later phases can reference
concrete names and real numbers instead of placeholders and assumptions.

## Phase 1 — Design the new protein sequence type

- Define a new protein sequence representation: a contiguous `uint8_t`
  buffer plus length, with a fixed alphabet-to-code mapping (documented in
  one place, e.g. a header/constants module) covering the 20 standard
  amino acids, gap/indel, and the ambiguity codes the current code already
  handles (check Phase 0 findings for which ones are supported today —
  don't drop support for symbols the old code accepts).
- Decide encode/decode functions: ASCII-in-FASTA -> `uint8_t` code on
  load, and code -> ASCII for any place that prints/serializes sequences
  (error messages, output files). These must be exhaustively covered by
  unit tests for every symbol in the alphabet, including lower/upper case
  and any "unknown residue" fallback behavior the old code has.
- Decide the runtime/compile-time dispatch strategy for SIMD: e.g. build
  multiple code paths (scalar, SSE4, AVX2) with a runtime CPU-feature
  check selecting the fastest available at startup, or a compile-time flag
  if FastPhylo's build/deployment model makes that acceptable. Confirm
  which is consistent with how the DNA fast-distance code already handles
  this, if it does anything similar — reuse that mechanism if it exists
  rather than inventing a second one.
- Write this up as a short design note (a page is enough) before
  implementing, covering: the struct/class layout, the alphabet mapping
  table, the dispatch mechanism, and exactly which existing functions will
  be reimplemented vs. left untouched in this pass.

## Phase 2 — Implement the new comparison primitives

- Implement byte-array mismatch counting ("how many positions are
  different") using SIMD compare + popcount, with a correct scalar
  fallback for architectures without the relevant instruction set and for
  sequence-length remainders that don't fill a full SIMD register.
- Implement the amino-acid replacement/substitution count matrix
  computation (tally of `counts[a][b]` over aligned positions) as
  described in the design decisions above. Include the scalar version
  first, profile before adding any vectorized/histogram variant.
- Both functions must handle: sequences of unequal length gracefully
  (reject, truncate, or error — match whatever the existing code does),
  ambiguity/gap codes according to whatever semantics the current
  implementation uses (e.g. does a gap count as a mismatch? does an
  ambiguous residue contribute partial counts anywhere?), and empty
  sequences.
- Keep these as free functions or a small class with a narrow, testable
  interface — avoid entangling them with I/O or the wider distance-matrix
  pipeline at this stage; that wiring happens in Phase 3.

## Phase 3 — Wire into fastprot's distance computation

- Replace the call sites identified in Phase 0 (the `calculate_*_dists`
  functions in `ProtDistCalc.cpp`) with calls to the new primitives,
  behind the switch described in design decision 5, so both old and new
  paths remain reachable (e.g. an environment variable, a compile flag,
  or a function parameter threaded through for testing purposes only).
- `fastprot_mpi` shares the same `ProtSeqUtils.cpp`/`ProtDistCalc.cpp`
  logic (via its own copies) — decide whether to update it in lockstep or
  explicitly leave it on the old path for this round, and say which.
- fastphylo-py is out of scope (see "Target and scope" above): do not
  touch `fastphylo-py/src/cpp/bindings.cpp` or
  `fastphylo-py/src/fastphylo/protein.py` in this phase.
- Do not remove the old code path yet.

## Phase 4 — Correctness verification

Build a differential test harness that runs the old and new
implementations against the same inputs and asserts identical (or
tolerance-bounded, if floats are involved and the tolerance is justified
in writing) results. Cover:

- Randomly generated protein sequence pairs across a range of lengths
  (very short, typical ~200-500 aa, and long, e.g. several thousand
  residues), with mismatch rates spanning low to high.
- Every symbol in the alphabet (all 20 amino acids, gap, and each
  supported ambiguity code) exercised at least once, including at
  sequence boundaries (first/last position) and in runs (e.g. all-gap
  regions).
- Edge cases: empty sequence, single-residue sequence, identical
  sequences, completely different sequences, sequences that are all one
  symbol, sequence-length mismatches (whatever the defined error/handling
  behavior is).
- Real biological data: `examples/protein_seq.fasta`/
  `protein_seq.phylip` are the existing FastPhylo fixtures (used by
  `RunExamples.sh` examples 6/7) — checked: only 4 sequences of ~12
  residues each, far too small for the length/dataset-size ranges Phase 5
  needs. Use them for the correctness smoke test (they're the existing
  convention), but pull in a small public MSA (real protein family,
  hundreds of residues, dozens+ of sequences) for the meaningful
  correctness and benchmark runs.
- The full end-to-end distance-matrix output (not just the low-level
  mismatch/matrix primitives) for at least one real dataset, comparing
  old-pipeline output to new-pipeline output — i.e. `fastprot`'s actual
  CLI output, in the spirit of `RunExamples.sh`'s diff-against-
  `expected_output/` convention.
- Wire this differential harness into FastPhylo's own test/example
  infrastructure so it runs repeatably: extend `RunExamples.sh` and/or
  get `code_tests/` actually building and running (per the "no CI today"
  finding in the scope note above) rather than assuming a CI pipeline to
  hook into. If there's no reasonable way to run it automatically yet,
  say so plainly rather than calling it "wired into CI."

This phase does not pass until every category above has an automated,
repeatable test that is checked into the repo.

## Phase 5 — Performance benchmark

Build a benchmark that is honest about what it measures and resistant to
common pitfalls (JIT/cache warm-up, measurement noise, unrepresentative
input sizes). Requirements:

- Benchmark both primitives in isolation (mismatch counting; replacement
  matrix tally) and the end-to-end protein distance computation, old vs.
  new, on identical inputs.
- Vary sequence length (short/typical/long, matching real protein sizes)
  and dataset size (number of sequences, since FastPhylo's real workload
  is all-pairs or many-pairs distance matrices, not a single pair).
- Run enough repetitions to report a distribution (median and
  interquartile range, or mean ± stddev), not a single timing — discard or
  flag runs affected by obvious system noise. Include a warm-up phase
  before timed runs.
- Run on the actual target hardware/architecture(s) FastPhylo is meant to
  support, including a case without AVX2/AVX-512 if the scalar fallback
  path matters for real users — the benchmark should demonstrate the
  fallback path isn't a regression, and the SIMD path's speedup where
  available.
- Report results in a small, checked-in table/script (e.g. a benchmark
  script plus a results markdown/CSV), not just console output, so the
  speedup claim is reproducible by re-running one command.
- Explicitly state the acceptance bar for "faster" (e.g. a target speedup
  factor decided with Lasse before this phase, since "faster" alone isn't
  a pass/fail criterion) and report whether it was met, including for the
  replacement-matrix tally specifically — it may not benefit as much from
  this change as the mismatch counter, and that should be reported
  honestly rather than obscured by an aggregate number.

## Phase 6 — Cleanup and finalize

Only after Phase 4 and Phase 5 both pass and results have been reviewed:

- Remove the old protein comparison code path and the switch/flag used to
  select between old and new during development.
- Update any developer-facing documentation (README, code comments,
  fastphylo-py docstrings) describing the protein sequence representation
  and comparison functions.
- Update the CHANGELOG (or equivalent) noting the performance change and,
  if relevant, any minimum-CPU-feature requirement introduced.

## Future phases (explicitly out of scope for this plan)

Noted here so they aren't lost, but none of these should be started
until this plan's Phase 6 is done and signed off. Each is a separate
bottleneck from the byte-encoding/SIMD work above and deserves its own
plan rather than being folded in opportunistically.

- **Newton-Raphson solver for the ML protein models.** `fastprot`'s
  WAG/JTT/LG/... models (`calculate_ml_dists` in `ProtDistCalc.cpp`) use
  Newton-Raphson (`MaximumLikelihood.cpp:likelihood_calc`), not Brent —
  Brent lives only in fastphylo-py's out-of-scope `bindings.cpp`. Each
  Newton iteration calls `likelihood_deriv()`, which computes a full
  20x20 matrix exponential via `Matrix::expm()` twice per step (current
  value plus a finite-difference probe), for up to 50 iterations. This
  is a different kind of bottleneck than this plan's — numerical linear
  algebra cost per pair, not residue-level comparison — and doesn't
  affect the Kimura path at all (Kimura is closed-form, no solver
  involved). Likely candidates for a future pass: caching/reusing
  `expm()` results, an analytic derivative instead of finite differences,
  or a cheaper matrix-exponential algorithm for the fixed 20x20 case.
  Profile this path specifically (independent of this plan's Phase 0
  profiling of the Kimura/count-matrix path) before committing to an
  approach.
- **Turn FastPhylo into a library fastphylo-py links against**, instead
  of fastphylo-py vendoring/reimplementing core logic (see "Target and
  scope" above). Would remove the current duplication between
  `fastprot`'s C++ and fastphylo-py's independent `bindings.cpp`
  reimplementation, and let fixes/speedups made in FastPhylo core
  actually reach fastphylo-py users.
- **Modernize `fastdist`'s DNA SIMD path.** Currently hardcoded SSE2
  only (`-msse2`, no AVX2/AVX-512), no runtime CPU-feature dispatch (see
  "Target and scope" above). Real-world speedup from the current SIMD
  code is unclear and worth measuring before investing further.

## Final verification checklist (before calling this done)

- [ ] All new alphabet/encoding edge cases have passing unit tests.
- [ ] Differential test suite (old vs. new) passes on synthetic and real
      data, checked into CI.
- [ ] End-to-end distance-matrix output matches old implementation on at
      least one real dataset.
- [ ] Scalar fallback path is tested and correct on a machine/build
      without the assumed SIMD extension.
- [ ] Benchmark results are checked in, reproducible via a single command,
      and show whether the agreed speedup target was met — separately for
      mismatch counting, replacement-matrix tally, and end-to-end protein
      distance.
- [ ] Old code path fully removed, no dead flags/switches left behind.
- [ ] Documentation and changelog updated.
