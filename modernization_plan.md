# Plan: Modernize FastPhylo to C++17, modernize tooling

## Goal

Bring FastPhylo's code style up to modern C++ (C++17 floor) and modernize
its build/dev tooling, without changing any program's observable
behavior. Same discipline as the speed2026a work this branches from:
verify byte-identical output (`examples/RunExamples.sh`, `ctest`) after
every phase, never trade correctness for style.

Branch: `modernize-cpp17` (off `speed2026a`'s tip, so it carries forward
all the performance work rather than starting from `master`).

## Scope for this plan

Confirmed with Lasse: **core library + `fastprot` first**, other
programs (`fastdist`, `fnj`, `fastprot_mpi`, `buildtree`,
`CreateSimulatedData`, `sequence_nj`) as later, separate phases. Six
standalone programs identified as no-longer-needed were removed
outright rather than modernized: `Tamir2Phylip.cpp`,
`Clustal2gaplessPhylip.cpp`, `drawtree.cpp`, `RF_dist.cpp`,
`naivedist.cpp`, `aml_tree.cpp` (done - see this branch's first two
commits).

**"Core library" scoping decision** (my call, flagging for review):
`libfastphylo.a` (`FASTPHYLO_SRCS` in `src/c++/CMakeLists.txt`) is
~14,000 lines across everything every program links against, but a
large chunk of that - `DNA_b128/`, `distance_methods/`,
`sequence_likelihood/` (~8,250 lines) - is DNA/tree-distance code
`fastprot` never touches at runtime, and nothing in this whole
engagement (SIMD work, output speedups, ML speedups) has read, tested,
or verified any of it. Proposing to **exclude those three directories
from this phase** and treat them as their own later phase, scoped and
audited on their own terms (`DNA_b128` in particular is SIMD-heavy code
that already got dedicated attention once, in the SSE2-on-ARM/simde fix
- it deserves the same care here, not a mechanical pass alongside
unrelated protein code). Say so if this should be wider.

**In scope for this phase** (~6,000-7,000 hand-written lines):
- `src/c++/`: `BitVector`, `Exception`, `InitAndPrintOn_utils`,
  `Object`, `Sequence`, `SequenceTree`(`_MostParsimonious`), `Tree`(`_impl`),
  `Simulator`, `arg_utils_ext`, `arg_utils` (C), `file_utils`,
  `stl_utils`, `std_c_utils` (C), `xml_output_global`, `DistanceMatrix`
  (`_impl`), `FloatDistanceMatrix` (`_impl`), `DistanceRow` (`_impl`),
  `log_utils.hpp`, `nucleotide.hpp`.
  (`fileFormatSchema.hpp` excluded - it's a generated/embedded XML
  schema string, not hand-written code to "modernize".)
- `src/c++/programs/fastprot/`: all hand-written sources (main.cpp,
  Prot*/ExpectedDistance/MaximumLikelihood/Matrix/*Stream* files).
  Gengetopt-generated code (`fastprot_gengetopt.*`, build-generated) is
  out of scope - it's not hand-written and regenerates from
  `gengetopt/cmdline.ggo` on every build.
- Build tooling for the whole project (CMake, clang-format/tidy, CI) -
  confirmed in scope, since it affects every program, not just this
  phase's file scope.

## Current-state audit

- **Language standard**: `src/c++/CMakeLists.txt` targets C++11 (and,
  until a bug fixed on `speed2026a`, wasn't even reliably reaching that
  on Clang - see that branch's first commit). No `CMAKE_CXX_STANDARD`
  property used anywhere; the flag is a hand-appended string.
- **CMake**: `cmake_minimum_required(VERSION 3.2)` (~2015-era). Old-style
  global commands throughout (`INCLUDE_DIRECTORIES`, `SET`, no
  `target_include_directories`/`target_compile_features` anywhere).
- **Tooling**: no `.clang-format`, no `.clang-tidy`, no CI
  (`.github/workflows` doesn't exist), no package manager config.
- **Code patterns** (rough greps across all of `src/c++`, not just this
  phase's scope - narrower counts TBD per-file during Phase 1):
  85 files with old-style `#ifndef`/`#define` include guards (0 use
  `#pragma once`), 61 with `NULL`, 44 with `typedef`, 24 with raw
  `new`/`delete`, a smaller number of C-style casts.
- **Legacy cruft found and removed this branch**: `Makefile.fromIsaac`,
  a dead pre-CMake build file from 2006 (kept, since `buildtree.cpp`/
  `fastdist.cpp`/`CreateSimulatedData.cpp` still nominally reference it,
  but its dangling references to the 6 removed programs were cleaned
  up).
- **Known correctness landmines already found in this codebase this
  engagement** (context for why every phase here gets verified, not
  assumed safe): a missing `#include` that left undefined behavior
  compiled into `fastprot` on optimized builds; `delete` called on a
  non-`new`'d pointer, corrupting the heap; the C++11 flag silently
  not applying on Clang. This is not a codebase where "obviously safe"
  mechanical changes can be trusted without building and testing.

## Phases

### Phase 0 - Tooling first [DONE, 2026-08-01]

Commits: "Phase 0: tooling..." and "Add basic CI...". Delivered:
`.clang-format` (LLVM-based, undocumented existing style so this is a
fresh baseline, not applied to existing files yet), `.clang-tidy`
(report-only, verified reasonable noise on both an old and a recent
file), `.gitignore` (none existed - `build/` had been cluttering every
`git status` all engagement), `CMAKE_CXX_STANDARD 17` replacing the old
hand-appended flag string, `cmake_minimum_required` bumped 3.2 -> 3.16,
and a first CI workflow (build + `ctest` + `RunExamples.sh`, Ubuntu,
unverified by an actual GitHub Actions run since that needs a push).

Getting to real C++17 (not just declared) required two **necessary**
fixes, both project-wide blockers rather than style choices: removing
`file_utils.hpp/.cpp`'s dynamic exception specifications
(`throw(Exception)`, removed entirely in C++17 - only occurrence
anywhere in the tree) and `DNA_b128/computeTAMURANEIDistance_
DNA_b128_String.cpp`'s `register` keyword (also removed in C++17) -
the latter despite `DNA_b128` being explicitly out of this phase's
style-modernization scope, because it's part of `libfastphylo.a`,
which every program links against, so it was blocking the whole build,
not just its own directory. Worth noting as a general lesson for later
phases: the C++17 *language standard* applies project-wide the moment
it's turned on, even though *style* modernization is scoped
incrementally - a file can be out of scope for style changes and still
need a minimal, required syntax fix to keep compiling.

Full rebuild + `RunExamples.sh` + `ctest` verified byte-identical
throughout (same discipline as speed2026a).

Not yet started: CMake target-based migration (`target_include_
directories`/`target_compile_features` for `fastphylo`/`fastprot`) -
Phase 0's plan included this but it wasn't reached this round; carrying
it forward rather than treating it as silently dropped.

### Phase 0 (original plan text, for reference)

Do this before any code changes, so subsequent phases get consistent
formatting/linting as they go instead of reformatting twice.

- Add `.clang-format` (a reasonable modern default - LLVM or Google
  base style, adjusted minimally; propose a starting point, iterate
  with Lasse rather than bikeshed alone).
- Add `.clang-tidy` with a conservative initial check set
  (modernize-*, readability-*, bugprone-*), run in report-only mode
  first to gauge noise before enabling anything as a build-breaking
  gate.
- CMake: bump `cmake_minimum_required`, add `CMAKE_CXX_STANDARD 17` /
  `CMAKE_CXX_STANDARD_REQUIRED ON` (replacing the hand-appended flag
  string), begin migrating the most heavily-touched targets
  (`fastphylo`, `fastprot`) to `target_include_directories`/
  `target_compile_features` - full CMake modernization of every target
  in the project is a larger, separate follow-on, not blocking this
  phase.
- Add a basic CI workflow (`.github/workflows/`): build + run
  `RunExamples.sh` + `ctest` on push/PR. There is none today - this is
  new safety net, not a replacement for local verification.

### Phase 1 - Mechanical, low-risk modernization

Changes with no semantic ambiguity, verifiable by build + full
regression pass, ideally partially automatable via `clang-tidy --fix`
once Phase 0's config exists (verify its suggested diffs before trusting
them, per this codebase's track record):

- `#ifndef`/`#define` include guards -> `#pragma once`.
- `NULL` -> `nullptr`.
- C-style casts -> `static_cast`/`reinterpret_cast`/`const_cast` as
  appropriate (verify each - a cast doing an implicit narrowing or
  pointer-type change needs the right one, not a mechanical swap).
- `typedef` -> `using`.
- `0`/`false` used where `nullptr` is meant (pointer comparisons/
  assignments specifically, not integer 0).

### Phase 2 - RAII / ownership

Higher risk, needs per-site review, not blanket automation:

- Raw `new`/`delete` (24 files) -> `std::unique_ptr`/`std::vector`/RAII
  wrappers, case by case. Verify ownership/lifetime semantics don't
  change (a `delete` in a destructor vs. a `delete` transferring
  ownership elsewhere are different fixes).
- Given this session already found one heap-corruption bug from
  exactly this class of code (`delete` on a non-`new`'d `FILE*`), audit
  every remaining raw `new`/`delete` site for similar mismatches while
  touching it, not just style.

### Phase 3 - Modern idioms

- Iterator-loop patterns (`for (it = x.begin(); it != x.end(); it++)`,
  common throughout this codebase per direct observation this session)
  -> range-based `for` where the loop body doesn't need the iterator
  itself.
- `override`/`= default`/`= delete` on virtual/special member functions
  where applicable.
- `auto` where it improves readability without obscuring an
  important type.
- Algorithm-header usage (`<algorithm>`) in place of manual loops where
  it's a clear readability win, not a stylistic detour.

### Phase 4 - Verification and sign-off

- Full `RunExamples.sh` + `ctest` pass, same as every phase.
- `benchmarks/run_benchmarks.sh` re-run to confirm modernization didn't
  regress speed2026a's performance work (style changes shouldn't, but
  verify rather than assume, especially anything touching `Matrix`/
  `ProtSeqCompare`/hot loops).
- Short writeup per phase (this file, updated) - same convention as
  `phase0_audit.md` for the performance work.

## What's explicitly deferred, not forgotten

- `DNA_b128`/`distance_methods`/`sequence_likelihood` (own future
  phase, per the scoping decision above).
- `fastdist`, `fnj`, `fastprot_mpi`, and the smaller remaining programs
  (`buildtree`, `CreateSimulatedData`, `sequence_nj`) - not touched this
  phase.
- The `Matrix` class's heap-allocation-per-instance design (flagged
  separately in memory - `MatrixExpm`/Eigen investigation) - that's a
  performance-motivated redesign, distinct from this plan's style
  modernization, though Phase 2's RAII work here may touch related code.
- Full CMake target-based migration for every target (`fastdist`, `fnj`,
  `fastprot_mpi`, docbook) - Phase 0 only starts this for `fastphylo`/
  `fastprot`.
