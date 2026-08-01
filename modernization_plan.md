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

### Phase 1 - Mechanical, low-risk modernization [DONE, 2026-08-01]

Commits: "Phase 1: convert include guards...", "Phase 1: NULL ->
nullptr...", "Phase 1: typedef -> using...", "Phase 1: C-style
casts...", and "Fix latent UB: Object::equals()...". Delivered:

- 40 headers converted from `#ifndef`/`#define` guards to
  `#pragma once`, via a small script (not blind regex) that tracks
  `#if`/`#ifdef`/`#ifndef` nesting depth to find the guard's *matching*
  `#endif` - verified against the most nested file (40 conditional
  directives) before running broadly. Also had to handle one
  non-UTF-8-encoded file (a legacy Latin-1 byte in a comment) via
  lossless round-tripping rather than erroring or corrupting it.
- 138 `NULL` -> `nullptr` (18 files) and 19 C-style casts -> the
  correct `static_cast`/`reinterpret_cast`/`const_cast`/`dynamic_cast`
  per clang-tidy's own categorization of each finding, not a blind
  `static_cast` sweep.
- 37 `typedef` -> `using` (14 files) - done **by hand**, not
  `clang-tidy --fix`: the automated fix-it corrupted 5 of the 37 sites,
  pasting in unrelated surrounding source (a class body, a commented-out
  code block) instead of the intended one-line replacement, when
  generating replacements near macros or unusual comments. Caught by
  the full rebuild (undefined-type errors), not by trusting the tool's
  output. Two anonymous typedef-struct/typedef-enum idioms were given
  real names directly (the more idiomatic modern-C++ equivalent) rather
  than mechanically wrapped in a `using` - deliberately did *not*
  convert the enum to a scoped `enum class`, since that changes every
  call site referencing its values, which belongs in Phase 3.
- **Found and fixed a real, disclosed bug along the way** (not silently
  bundled into the "no behavior change" cast-modernization commit):
  `BitVector::equals()`/`Sequence::equals()` downcast `const Object*`
  with a C-style cast, which for a base-to-derived polymorphic
  conversion is an *unchecked* cast - comparing to a differently-typed
  `Object` was undefined behavior, not a defined "not equal". Fixed
  with `dynamic_cast` (BitVector's existing null check was dead code
  before this; added the equivalent check to Sequence, which didn't
  have one at all). Verified with a standalone cross-type comparison
  test, not just "the existing tests still pass" (they don't exercise
  this path at all).

Lesson reinforced for later phases: don't trust `clang-tidy --fix`
output blindly on this codebase, even for checks that look purely
mechanical - verify by rebuilding, not by reading the diff alone.

### Phase 2 - RAII / ownership [DONE, 2026-08-01]

Commits: "Phase 2 (RAII): FloatDistanceMatrix owns rows by value...",
"Phase 2 (RAII): BinaryDmOutputStream owns its file stream...", "Phase 2
(RAII): fastprot main.cpp owns its input/output streams...". Delivered,
all three raw `new`/`delete` sites found in this phase's file scope:

- **`FloatDistanceMatrix`**: each row was `std::vector<DistanceType>*`
  (manually `new`'d in `assureSize()`, `delete`'d in the destructor and
  `removeLastRow()`), with no actual need for the indirection - nothing
  aliased a row's address, no polymorphism, no optional rows. Changed to
  a plain `std::vector<std::vector<DistanceType>>`. This also fixed a
  **real, previously-undiscovered bug**: `operator=()` resized `D` and
  immediately dereferenced the new (nullptr-initialized) row pointers
  without allocating them first, so copy-constructing or copy-assigning
  a `FloatDistanceMatrix` crashed unconditionally on a null-pointer
  dereference. The one call site
  (`distance_methods/LeastSquaresFit.cpp:186`) turned out to be
  unreachable from any built program, which is presumably why nobody had
  hit it - but it's a live bug in the library's public API, not
  theoretical. Removed the hand-written copy constructor/`operator=`
  entirely (Rule of Zero - correct by construction once every member
  manages its own storage) rather than just patching the null-deref.
  Verified with a standalone test exercising copy-construction/
  assignment directly, since neither `RunExamples.sh` nor `ctest`
  reaches this path.
- **`BinaryDmOutputStream`** (fastprot): `ofs` was a raw `ostream*` that
  sometimes aliased `&cout` (not owned) and sometimes pointed to a
  heap-allocated `ofstream` (owned, tracked with a separate `bool` and
  manually `delete`'d). Split ownership from aliasing instead: a new
  `std::unique_ptr<std::ofstream> file_stream` member owns the file, and
  `ofs` stays a non-owning observer pointer. Hand-written destructor
  removed entirely. `open_write_binary()` (`file_utils.cpp`) still
  returns a raw `ofstream*` - left as-is since `fastdist`'s and
  `fastprot_mpi`'s `BinaryDmOutputStream` (out of this plan's scope)
  also call it and assign straight into a raw `ostream*`; changing its
  signature would force out-of-scope files to change too.
- **`fastprot/main.cpp`**: `DataInputStream*`/`DataOutputStream*` were
  `new`'d in a format-selection `switch` and `delete`'d at the end of
  the `try` block - except the `catch(...){ throw; }` handler re-throws
  unconditionally, so any exception during read/write (a malformed
  input file, a full disk) skipped the `delete` calls. Real leak on
  every error path, not just a theoretical one, given this code parses
  untrusted user-supplied files. Switched both to `std::unique_ptr`
  (`std::make_unique` at each assignment site); cleanup is now automatic
  on every path, including the ones that used to leak.

Same verification discipline throughout: full rebuild + `RunExamples.sh`
+ `ctest` byte-identical after each commit (two pre-existing, unrelated
`ex4`/`ex10` mismatches confirmed present before this phase too, not
caused by it).

### Phase 3 - Modern idioms [DONE, 2026-08-01]

Commits: "Phase 3: convert manual iterator for-loops to range-based
for...", "Phase 3: annotate virtual overrides with override...",
"Phase 3: use auto where it avoids duplicating a static_cast's target
type...". Delivered:

- **Range-based for**: 6 sites (`SequenceTree.cpp` x3, `stl_utils.cpp`,
  `ExpectedDistance.cpp` x2, `Matrix.cpp`, `ProtSeqUtils.cpp` x2) where
  the classic `for (it = x.begin(); it != x.end(); it++)` loop body
  only ever dereferenced the iterator. Also removed the now-unnecessary
  iterator variables that had been hoisted above their loop. Left one
  loop (`ProtSeqUtils.cpp`'s gap-removal, iterating `positions.rbegin()`
  /`rend()`) as a manual loop with `auto` for the iterator type instead
  of range-based for, since it deliberately iterates in reverse (erasing
  higher string positions before lower ones so indices don't shift) -
  there's no range-based for syntax for reverse iteration without a
  wrapper type.
- **`override`**: added to every virtual function in scope that
  genuinely overrides a base virtual - `Object`'s
  printOn/objInitFromStream/hashCode/equals/destructor as overridden by
  `BitVector`, `Sequence`, `Exception`, `DistanceMatrix`,
  `FloatDistanceMatrix`, `Tree`, `TreeNode`; and
  `DataInputStream`'s/`DataOutputStream`'s
  destructor/read/readSequences/print/printSD/printStartRun/printEndRun
  /printHeader as overridden by the six fastprot stream classes.
  Deliberately did **not** add it to virtual functions a derived class
  merely *introduces* (e.g. `Sequence::printShort`/`printWithoutGaps`
  have no base declaration to override) or to the base classes' own
  first declarations - getting this right required reading each base
  class, not pattern-matching "virtual" near "public Base". Full
  rebuild (including `fastdist`/`fnj`, which share the touched
  core-library headers but are out of this plan's scope) came back
  clean with no `-Winconsistent-missing-override` warnings.
- **`auto`**: applied only where it removes a duplicated type name (6
  sites, `unsigned char up = static_cast<unsigned char>(...)` ->
  `auto up = static_cast<unsigned char>(...)`, in speed2026a's SIMD
  comparison code). Deliberately left `Sequence.cpp`'s
  `dynamic_cast`-typed variable alone (the explicit type is load-bearing
  context for the comment explaining that bug fix) and skipped a
  `double matches = static_cast<double>(a) - static_cast<double>(b);`
  site (a computed expression, not a bare same-type cast, so `auto`
  wouldn't remove any duplication) - matches the plan's caution about
  not obscuring an important type.
- **`<algorithm>`**: searched the full phase scope; found no manual
  loops that were a clear readability win to replace beyond what
  speed2026a's work already introduced (`std::min`,
  `std::accumulate`). Not forcing conversions here, per the plan's own
  "not a stylistic detour" caution.

Same verification discipline throughout: full rebuild + `RunExamples.sh`
+ `ctest` byte-identical after each commit (same two pre-existing,
unrelated `ex4`/`ex10` mismatches as every prior phase - see Phase 4,
which finally ran these down: not a real bug).

### Phase 4 - Verification and sign-off [DONE, 2026-08-01]

- **Full rebuild from scratch**: deleted the build directory entirely
  and reconfigured/rebuilt in `Release` mode (every prior phase's
  verification reused an incrementally-updated `Debug` build) - clean,
  no errors or warnings, confirming nothing in Phases 0-3 was
  accidentally relying on stale build state. `ctest` and
  `RunExamples.sh` both pass against this from-scratch build,
  byte-identical to `expected_output/`.
- **Finally ran down the "ex4/ex10 mismatch" this whole branch kept
  reporting as pre-existing-but-unexplained**: they're not a real
  output mismatch at all. Both `run_example` calls are commented out in
  `examples/RunExamples.sh` (`#run_example 4 ...`, `#run_example 10
  ...`) - predates this branch, unrelated to any of this work - so no
  `ex4.out`/`ex10.out` is ever generated. Every verification pass this
  branch compared `expected_output/ex{4,10}.out` (which still exist as
  fixture files) against a nonexistent generated file, which trivially
  "differs" regardless of code correctness. The script's own built-in
  check (`for i in ex*.out; do diff ...`) is immune to this since it
  only globs files that were actually generated - the bug was in this
  session's own ad hoc verification loop, not in the codebase.
  Corrected description throughout this doc; the accurate statement for
  every phase this branch is: **all 15 active examples byte-identical,
  always** - there was never an unexplained mismatch to carry forward.
- **`benchmarks/run_benchmarks.sh` re-run** to confirm the style
  modernization didn't regress speed2026a's SIMD work. Primitive-level
  absolute numbers (`bench_primitives`, which still has its own
  internal scalar-vs-SIMD comparison, unrelated to fastprot's) came out
  slower in wall-clock terms than the checked-in Phase 5 baseline
  (e.g. `count_id_dist` at length 300: 2.375µs here vs 1.625µs
  checked-in) - but the *old* column scaled up by essentially the same
  factor (16.8µs vs 11.5µs checked-in), and the resulting speedup
  ratios matched the checked-in baseline closely (7.07x vs 7.10x, etc.)
  - consistent with general machine load during this run, not a code
  regression, since a real regression would only slow the *new*
  column. End-to-end `fastprot -D JCK` wall-clock times (the number
  that actually matters for real usage) came out on par with or
  slightly *faster* than the Phase 5 baseline at all four problem
  sizes. Note: the old-vs-new toggle inside `fastprot` itself
  (`FASTPROT_NEW_PROT_CODE`) no longer exists - it was removed in
  speed2026a's own Phase 6 cleanup - so the end-to-end script's "old"
  and "new" columns are now identical code measured twice; only the
  "new" column (or either, they're the same) is meaningful post-Phase-6.
  Did not overwrite the checked-in `results_primitives.csv`/
  `results_end_to_end.csv` with this run's noisier numbers - reverted
  them after comparing, since the checked-in files are the
  well-documented, controlled Phase 5 measurement and this was a
  confirmatory spot-check, not a new canonical baseline.

**Overall sign-off for this plan's four phases**: Phase 0 (tooling),
Phase 1 (mechanical: pragma once/nullptr/using/casts), Phase 2 (RAII:
3 raw-ownership sites converted, 2 real latent bugs found and fixed
along the way - the `Object::equals()` UB and
`FloatDistanceMatrix::operator=()`'s null-deref), and Phase 3 (modern
idioms: range-based for/override/auto) are all complete for this
plan's scope (core library minus `DNA_b128`/`distance_methods`/
`sequence_likelihood`, plus `fastprot`). No behavior change anywhere
except the two disclosed, intentional bug fixes; every commit verified
independently. What's next is explicitly out of this plan's scope, see
below.

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
