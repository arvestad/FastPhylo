# Plan: fastphylo-py consumes libfastphylo instead of vendoring/reimplementing it

## Context

The `eigen-ml-matrix-backend` branch's protein maximum-likelihood (ML)
distance path was rebuilt specifically with a future Python-bindings
caller in mind: `MatrixExpm build_ml_decomposition(model_type mt)` does
the expensive one-time eigendecomposition step, and `calculate_ml_dists
(sv, dm, const MatrixExpm&)` takes an already-built decomposition and
has **no** overload that builds its own - a caller with an unknown/
unbounded number of distance-estimation calls (exactly what a
long-lived Python session is) can never accidentally redo the
expensive step. This plan makes fastphylo-py (github.com/arvestad/
fastphylo-py, PyPI package `fastphylo`, pybind11 + scikit-build-core,
GPLv3, MIT-compatible) an actual consumer of that work, instead of its
current two independent problems:

1. `src/cpp/fastphylo/` (52 files) is a **stale, vendored fork** of
   this repo's pre-restructuring core+DNA code - confirmed still has
   `Object.hpp/.cpp` and `Exception.hpp/.cpp`, both deleted from real
   FastPhylo during this engagement's legacy-error-handling cleanup.
2. `src/cpp/bindings.cpp` has an **independent, from-scratch
   reimplementation** of protein ML distance estimation (its own
   count-matrix builder, its own Brent's-method optimizer) that never
   calls FastPhylo's code at all. Its Python side (`protein.py`)
   already does its own eigendecomposition via `numpy.linalg.eigh` on
   the *same* symmetrized similarity-transform trick this session
   built into `MatrixExpm` in C++ - independent validation the
   approach is sound.

Goal: make fastphylo-py a real consumer of FastPhylo's C++ core for
both the protein ML path and the DNA/tree path, preserving its
existing Python-facing API (`distance_matrix()`, `nj()`/`fnj()`/
`bionj()`, `Tree`, `DistanceMatrix`) unchanged.

## Decided (not open)

- Repos stay separate; FastPhylo is pulled into fastphylo-py via a git
  submodule pinned to a tagged FastPhylo release commit.
- Every protein model fastphylo-py supports moves into FastPhylo's C++
  `ModelMatrix.cpp`, not just the 5 it had before this plan.
- A pure-Python fallback path is kept for protein ML when the compiled
  `_fastphylo` extension isn't available.
- Speed verification is a required, dedicated deliverable (Python-API-
  level measurement, not extrapolated from FastPhylo's own C++
  benchmarks), covered in a later phase once fastphylo-py's own
  bindings exist to measure.

**Explicitly not part of this plan** (tracked separately): `ExpectedDistance.hpp`'s
file-scope mutable static state (`log_prior_dist`, `prior_prob`,
`nr_distances`, `std_deviation`, etc.) is a real long-lived-caller/
thread-safety concern, but fastphylo-py's expected-distance path is
pure Python today and never touches this C++ code - nothing here
depends on fixing it.

## Phase 0 - Port fastphylo-py's 10 non-native protein models into `ModelMatrix` (DONE 2026-08-11)

Confirmed directly from fastphylo-py's real `src/fastphylo/protein.py`:
14 models total, not 16 as an earlier, less careful lookup suggested.
Overlap with FastPhylo's native 5 (`wag, day, jtt, mvr, lg`) is 4, not
5 - `mvr` (Müller-Vingron) has no fastphylo-py counterpart at all. The
real port list was **10 models**: `JTT-DCMUT, VT, HIVb, HIVw, cpREV,
BLOSUM62, DCMUT, MtREV, RtREV, PMB`.

fastphylo-py's amino-acid ordering (`ARNDCQEGHILKMFPSTWYV`) is
character-for-character identical to FastPhylo C++'s - no index-
permutation step needed when transcribing matrices.

**Data source and verification**: primary source of truth was
fastphylo-py's own `R` (190-value upper-triangular exchangeability) +
`freq` (20-value equilibrium) literals per model, cross-checked against
IQ-TREE2's `model/modelprotein.cpp` (same procedure that caught the
JTT 100x bug and the LG scale bug earlier this branch). Result: **9 of
10 models matched exactly** (DCMUT's R values are IQ-TREE2's data
times 0.01 - a harmless scale-convention difference, same class as
LG's, handled automatically by `mean_substitution_rate()`-based
normalization, not a bug). **VT is the one exception**: fastphylo-py's
VT data does not match IQ-TREE2's published VT numbers - not a scale
factor, both R and equilibrium frequencies differ by a few percent.
Confirmed VT's own data is still a valid, self-consistent rate matrix
(row sums to zero at machine precision, `1e-17` order), so this isn't
a transcription bug like JTT's was - most likely two different
published VT parameterizations in circulation. fastphylo-py's own
numbers were used (not IQ-TREE2's) to preserve exact behavioral
continuity for existing fastphylo-py users. Disclosed here rather than
silently picking one source.

**Storage format**: rather than transcribing a literal 400-entry `Q`
per model the way `get_wag()` etc. do, the 10 new models are stored as
`R` (190 values) + `eq` (20 values), with `Q(i,j) = R(i,j)*eq[j]`
off-diagonal and `Q(i,i) = -row_sum` computed, not hand-typed - this
directly targets the JTT bug's actual failure mode (a hand-typed
diagonal value can be silently wrong; a computed one cannot).
`build_ml_decomposition()` needed no changes at all - it's already
scale-agnostic per model, exactly the LG fix's mechanism.

**Extensibility**: `get_model_matrix()`/`get_model_vec()` now dispatch
the 10 new models through a `find_model_entry()` table lookup
(`{model_type, name, R*, eq*}` entries) rather than N separate
functions - the original 5 models' existing per-model getter functions
were left untouched (lower risk than retrofitting proven-correct data
to a new format). Adding model #16 later means one table entry, not a
new function pair plus new switch cases in two places.

Wired into `tests/ModelMatrix_test.cpp` (all 15 models, row-sum/
frequency-sum invariants) and `tests/MaximumLikelihood_test.cpp` (all
15 models' solver self-consistency against `examples/globin_family.fasta`)
in the same change, plus `apps/fastprot/main.cpp`'s `-D` flag and
`ProtDistCalc.cpp`'s dispatch switches (a secondary, low-cost benefit -
`fastprot` users get 10 new models too).

Verified: `ctest` 6/6 on Clang Release and real GCC 14;
`RunExamples.sh`/`RunCliChecks.sh` clean; all 10 new models produce
sensible ML output via the real `fastprot` CLI (HIVb/MtREV correctly
report "too diverged" on the globin-family fixture - expected, since
those are specialized models trained on very different sequence types).

## Phase 1 - FastPhylo: make `fastphylo` a real linkable library (DONE 2026-08-11)

1. Fixed the undeclared BLAS/LAPACK/Threads/libm dependency -
   `find_package`d and linked `PUBLIC` onto the `fastphylo` target
   itself (previously only attached to the executable targets, even
   though `Matrix.cpp` - part of `fastphylo`'s own sources - calls
   LAPACK directly). A consumer linking only `fastphylo` used to hit
   undefined-symbol errors; now doesn't.
2. Added `POSITION_INDEPENDENT_CODE ON` on the `fastphylo` target - a
   pybind11 extension module is itself a shared object; a non-PIC
   static archive can't be linked into one on Linux.
3. Added `FASTPHYLO_BUILD_APPS`/`FASTPHYLO_BUILD_TESTS` options
   (default `ON`, preserving native build/release behavior unchanged) -
   gate the 3 CLI apps and the 5 test executables respectively, so an
   external consumer embedding just the library doesn't waste build
   time compiling either.
4. Guarded the root `CMakeLists.txt`'s `enable_testing()`/`add_subdirectory
   (docs)`/CPack block behind `CMAKE_SOURCE_DIR STREQUAL
   CMAKE_CURRENT_SOURCE_DIR` (this repo's CMake floor is 3.18, no
   `PROJECT_IS_TOP_LEVEL` until 3.21) - these only make sense when
   FastPhylo is the top-level project, not pulled in as a subdirectory.
5. **Found and fixed a real bug this same work surfaced**: every path-
   construction variable in `src/c++/CMakeLists.txt`
   (`FASTPHYLO_CORE_DIR` etc., and several `target_include_directories`/
   `target_compile_definitions` calls) used `CMAKE_SOURCE_DIR`, which
   always resolves to the *outermost* project's source directory, not
   FastPhylo's own - harmless as long as FastPhylo was always the
   top-level project (true for every build this repo has ever done
   until this plan), but breaks immediately and completely when
   embedded via `add_subdirectory()` from an external project (every
   source file path resolves relative to the *external* project's
   directory instead). Fixed by replacing every `CMAKE_SOURCE_DIR`
   with `PROJECT_SOURCE_DIR` (correctly follows the nearest enclosing
   `project()` call, i.e. FastPhylo's own, regardless of nesting) in
   both `CMakeLists.txt` and `src/c++/CMakeLists.txt`. Caught by
   actually building a throwaway external "smoke test" project rather
   than assuming the CMake fixes above were sufficient - they weren't.
   Also fixed `vcpkg.json`'s missing `libxml2` dependency while in that
   area (not required for fastphylo-py, which does its own FASTA/
   Phylip/Stockholm parsing, but a pre-existing gap worth closing since
   the file was touched anyway).

**Checkpoint**: existing 3-platform-equivalent native build stayed
green throughout (`ctest` 6/6 on Clang Release and real GCC 14,
`RunExamples.sh`/`RunCliChecks.sh` clean) after every change. A
throwaway smoke-test CMake project (not committed - scratch only)
`add_subdirectory()`d this repo with both new options `OFF`, linked a
trivial executable against only `fastphylo`, and called one function
from each of core (`Sequence`/`StrDblMatrix`), dna (`computeNJTree`),
and protein (`build_ml_decomposition`/`calculate_ml_dists`) - confirmed
working end to end only after the `PROJECT_SOURCE_DIR` fix above.

### Side finding: a real, pre-existing crash bug in NJ/BioNJ for N=2 taxa (2026-08-11, not fixed here)

While building and iterating on the smoke test above (initially using
a 2-taxon input, the simplest possible case), `computeNJTree(dm, tree,
NJ)` crashed intermittently (segfault or bus error, ~90% of runs) -
confirmed **100% reproducible and unrelated to this plan's own changes**:
reproduces identically via the real `fnj -I phylip -m NJ` CLI on any
real 2-taxon PHYLIP file, and via a version of the same test compiled
directly against the pre-existing source tree, bypassing every CMake
change in this plan entirely. Root cause (via `lldb`, `EXC_BAD_ACCESS`
at `NeighborJoining.hpp:175`): `computeNeighborJoiningTree()`'s
`while (numNodes > 3)` loop is skipped entirely when there are only 2
starting taxa, falling straight through to an unconditional "resolve
the final 3-taxon star" block that reads `dm.getIdentifier(2)`/
`dm.getDistance(_, 2)` - out of bounds for a 2-row matrix.
`computeBioNJTree()` has the identical pattern (same fix needed).
`computeFNJTree()` is unaffected - by design it never computes final
branch lengths at all (`fnj`'s own default method, which is *why* this
was never noticed via the CLI's default invocation: `fnj`'s default
`--method` is `FNJ`, not `NJ`).

Not fixed as part of this plan - genuinely pre-existing, orthogonal to
the CMake/model-porting work here, and deserves its own careful fix +
verification rather than a rushed patch. The smoke test itself was
changed to use 3 taxa (also a more realistic scenario) to complete
Phase 1's checkpoint. Flagged here so it isn't lost; fix separately
when prioritized.

## Remaining phases (not started)

See the session record for the full remaining plan (Phase 2: consumption
mechanism/submodule wiring on the fastphylo-py side; Phase 3: pybind11
bindings redesign; Phase 4: pure-Python fallback; Phase 5: sequencing;
Phase 6: verification including the required speed benchmarks) - these
all require access to the fastphylo-py repository, which this session
did not have a local checkout of.
