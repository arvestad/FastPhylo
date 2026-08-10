# Changelog

## 2.0.0-beta.1 (Unreleased)

The first release since 1.0.3. Everything below accumulated on the
`modernize-cpp17` branch: a full C++17/CMake modernization, a project
layout restructure, a CLI11-based command line replacing gengetopt, a
cross-platform (Linux/macOS/Windows) CI and release pipeline, the
`speed2026a` SIMD performance work, and a long list of real bugs found
and fixed along the way. Given how much changed under the hood, this
is a beta - please report anything that looks wrong.

### Added

- Cross-platform (Linux/macOS/Windows) CI, and a release pipeline that
  builds, tests, packages, and publishes binaries for all three
  platforms from a pushed version tag (see `RELEASING.md`).
- `fnj` can now read the binary distance-matrix format it already
  wrote (previously write-only).
- A repo-wide code style (4-space indent, Allman braces, enforced by
  `clang-format` both locally via a pre-commit hook and in CI).

### Changed

- **Command line**: all four programs (`fastdist`, `fnj`, `fastprot`,
  `fastprot_mpi`) now parse arguments with CLI11 instead of gengetopt.
  Functionally equivalent flags throughout, but `--help`/`--version`
  output and error messages for invalid input are reformatted.
- **Build requirements**: CMake 3.18+ (was 3.2) and a C++17 compiler
  are now required project-wide.
- **libxml2 is now optional and auto-detected**, rather than a hard
  requirement. If it isn't found, `cmake` prints a warning and builds
  without XML support instead of failing; pass `-DWITH_LIBXML=OFF` to
  build without it intentionally and silence the warning.
  Correspondingly, `fastprot`/`fastdist -O` and `fnj -I` now default
  to Phylip format rather than XML.
- The old self-built `STATIC` toolchain (gengetopt/libxml2/zlib built
  from 2009-2011-era FTP-hosted tarballs) is gone; dependency
  detection now goes through CMake's standard `find_package`/vcpkg
  mechanisms.
- Project layout restructured around a conventional library shape:
  `include/fastphylo` (public API) + `src/fastphylo` (implementation,
  mirroring `include/fastphylo`) + `src/c++/apps` (the four
  command-line programs, each a thin consumer of the library).

### Performance

- `fastprot`'s ID/JC/JCK (Kimura)/JCSS distance models now use a
  byte-encoded, SIMD-accelerated comparison primitive
  (`ProtSeqCode`/`ProtSeqCompare`) instead of a per-character
  `toupper()`+compare loop. Measured 5-7x speedup on the comparison
  primitive itself and 1.8-6x end-to-end depending on dataset size
  (see `benchmarks/RESULTS.md` for full methodology and numbers,
  `planning/phase0_audit.md`/`planning/phase1_design.md` for the profiling and design
  work behind it). No minimum CPU feature is required beyond each
  platform's SIMD baseline (SSE2 on x86-64, NEON on AArch64, both
  present unconditionally on those architectures); a scalar fallback
  covers anything else, selected at compile time - no runtime
  dispatch. `count_replacements()` (used by the WAG/JTT/DAY/ARVE/MVR/LG
  expected-distance and maximum-likelihood models) was not in scope
  for this round.
- Matrix output (Phylip, XML, and binary) now writes one row per I/O
  call instead of one call per matrix entry - roughly 3-4% faster
  (Phylip), 2.5-3% (XML), up to ~17% on a print-heavy dataset.
- The maximum-likelihood protein distance code caches its rate
  matrix's eigendecomposition instead of recomputing it on every
  Newton-Raphson iteration, and cut temporary `Matrix` allocations
  from 6 to 3 per call.

### Fixed

- `fastdist -D JC`/`-D K2P`/`-D TN93` used to silently produce `NaN`/
  `Inf` for maximally-diverged sequence pairs (once the raw Hamming
  distance passes each model's mathematical domain limit - p≥0.75 for
  Jukes-Cantor, and analogous limits for Kimura/Tamura-Nei), which a
  downstream check then replaced with a generic `-1` and a "float not
  finite" warning - no indication of which pair or why. Each
  correction formula (including F84, not exposed via `-D` today but
  sharing the same code) now detects its own domain limit directly and
  reports a specific warning ("sequences too different for `<model>`
  distance, using distance=2.5") with a clamped, clearly-flagged
  distance instead.
- All four programs used to exit with status 0 even after printing a
  fatal error to stderr - any script checking the exit code was
  silently lied to.
- `fastprot`'s (and `fastprot_mpi`'s) `-O binary -o <file>` crashed on
  every invocation: freeing a `FILE*` that was opened with `fopen()`,
  not `new`. Binary output previously had no regression coverage at
  all; added.
- `fastprot`/`fastprot_mpi` crashed or misparsed arguments in
  optimized (`-O2`/`-O3`) builds on at least one platform (AppleClang
  on Apple Silicon): a missing `#include` left `WITH_LIBXML` undefined
  for the preprocessor, leaving a dead code path that read an
  uninitialized struct field - undefined behavior the optimizer
  exploited.
- The binary distance-matrix format had no run-boundary framing,
  silently corrupting output whenever `--number-of-runs`/`-r` > 1 was
  combined with `-O binary`. The wire format now records a run count
  so readers can tell separate runs apart.
- `fastdist`'s memory-efficient row-streaming path (`-e`) under
  `-D JC` printed zeros for the entire lower triangle instead of real
  distances (an inverted boolean check, out of step with its three
  sibling models).
- `fastdist` crashed on Phylip input with a trailing newline.
- Fixed an int/`size_t` mismatch in `DistanceMatrix` causing incorrect
  indexing on very large inputs, and a live memory-efficient-mode
  Jukes-Cantor bug, both found while consolidating 8 near-duplicate
  internal matrix-filling functions down to 3.
- `SequenceTree`'s Phylip reader enforced a stricter-than-necessary
  10-character name-field convention inherited from an older Phylip
  dialect, unable to parse some of this project's own plain fixture
  files.
- Fixed a self-assignment bug (`tree = tree` corrupting the tree via a
  use-after-free) and a variable-shadowing null-pointer bug in
  `makeCanonical` (certain single-node trees), plus an internal
  tree-validation check that ran an O(n) walk unconditionally, even in
  release builds.
- Numerous Windows/MSVC portability fixes (POSIX-only calls, VLAs, C99
  compound literals, dead inline assembly, stdout CRLF-translation
  corrupting binary output) - fastphylo now builds and passes its full
  regression suite on Windows, in addition to Linux and macOS.

### Removed

- The old hand-rolled `Exception` class, replaced by `std::runtime_error`.
- gengetopt and its self-built toolchain dependency.
- `Object.hpp`, a Java-style common base class every core data type
  used to inherit from for printing/equality/hashing. Replaced with
  idiomatic `operator<<`/`operator==`/`std::hash` - which also removed
  a real latent undefined-behavior risk (an unchecked pointer downcast
  in tree-equality comparison).
- Roughly 1200 lines of unused single-precision distance-matrix code
  (`FloatDistanceMatrix`), unused sequence-simulation code, several
  unused legacy standalone programs, and ~15 unused `SequenceTree`
  methods (including an unused Robinson-Foulds implementation) - none
  reachable from any of the four shipped programs.

### Known limitations

- `fastprot_mpi` (the MPI-distributed build of `fastprot`) is not
  covered by CI and has not been independently verified against the
  changes in this release. Building it remains best-effort only (see
  `INSTALL`); active development on it is paused for now.
