# Changelog

## Unreleased (speed2026a)

### Performance

- `fastprot`'s ID/JC/JCK (Kimura)/JCSS distance models now use a
  byte-encoded, SIMD-accelerated comparison primitive
  (`ProtSeqCode`/`ProtSeqCompare`, `src/c++/programs/fastprot/`) instead
  of a per-character `toupper()`+compare loop. Measured 5-7x speedup on
  the comparison primitive itself and 1.8-6x end-to-end depending on
  dataset size (see `benchmarks/RESULTS.md` for full methodology and
  numbers, `phase0_audit.md`/`phase1_design.md` for the profiling and
  design work behind it). No minimum CPU feature is required beyond each
  platform's SIMD baseline (SSE2 on x86-64, NEON on AArch64, both
  present unconditionally on those architectures); a scalar fallback
  covers anything else, selected at compile time - no runtime dispatch.
- `count_replacements()` (used by the WAG/JTT/DAY/ARVE/MVR/LG
  expected-distance and maximum-likelihood models) is unaffected - not
  in scope for this round (see `phase0_audit.md`'s scope notes).
- Output writing (`PhylipDmOutputStream::printPHYLIPfastSD`) is now the
  dominant cost for large distance matrices, having previously been
  masked by the slower comparison primitive - logged as follow-up work,
  not addressed here (see `phase0_audit.md`'s "Post-Phase-3 finding").

### Fixed

- `fastprot`/`fastprot_mpi` crashed or misparsed their arguments in
  optimized (`-O2`/`-O3`) builds on at least one platform (AppleClang on
  Apple Silicon): `main.cpp` was missing `#include "config.h"`, so
  `WITH_LIBXML` was never defined for the preprocessor there, leaving a
  dead code path that read an uninitialized struct field - undefined
  behavior that the optimizer exploited. Fixed by adding the include.
- `CMakeLists.txt` only applied `--std=c++11` `IF(CMAKE_COMPILER_IS_GNUCXX)`,
  never for Clang, so the project silently built in each compiler's
  pre-C++11 default there instead. Fixed by also matching Clang.
