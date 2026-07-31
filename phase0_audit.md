# Phase 0 audit — fastprot Kimura/protein-distance speedup

Scope: confirm `plan.md`'s claims against the real code, get a working
optimized build, and profile a realistic `fastprot` run. This note is the
Phase 0 deliverable.

## Blocking bug found and fixed first

`fastprot` could not be profiled at all in any optimized build
(`RelWithDebInfo`/`Release`) on this machine (Apple Silicon, AppleClang 15):
every invocation either printed a bogus "No input data or parameters" and
exited, or crashed with a SIGTRAP, regardless of the arguments given. A
plain `-O0` debug build worked fine.

Root cause: `src/c++/programs/fastprot/main.cpp` (and
`fastprot_mpi/main.cpp`) never `#include "config.h"` — unlike
`fastdist/main.cpp` and `fnj/main.cpp`, which do. `config.h` is where
`WITH_LIBXML` gets `#define`d (via `config.h.cmake` +
`CONFIGURE_FILE`); `WITH_LIBXML` is never passed as a compiler `-D` flag,
only used as a CMake-level variable. Without the include, every
`#ifdef WITH_LIBXML` / `#ifndef WITH_LIBXML` in `main.cpp` evaluates as if
libxml support were absent, even though the binary is in fact built with
it (`WITH_LIBXML:BOOL=ON`). Two consequences:

1. XML input support (`-I xml`) is silently dead code in the actual
   binary — the `#ifdef WITH_LIBXML` guarded `case input_format_arg_xml:`
   branch never compiles in.
2. The `#ifndef WITH_LIBXML` block *does* compile in, and reads
   `args_info.input_format_arg` — **before** `cmdline_parser()` has
   populated `args_info`, i.e. reading an uninitialized struct field.
   This is genuine UB. At `-O0` it happens to be harmless; at
   `-O2`/`-O3` Clang's optimizer exploits it, which manifested as the
   compiler folding the unrelated `isatty(STDIN_FILENO) && argc==1` check
   a few lines above into an unconditional branch (confirmed via
   disassembly — the compiled code calls `isatty` and then unconditionally
   falls into the "no input" exit path, with no compare/branch on `argc`
   at all) or, under slightly different codegen, a hard trap.

Fix: added `#include "config.h"` to both `main.cpp` files (one-line change
each, matching what `fastdist`/`fnj` already do). Verified:
- `examples/RunExamples.sh` (`FASTPHYLOPATH=build/src/c++`) — all
  examples, including `fastprot` examples 6/7, produce byte-identical
  output to `examples/expected_output/` before and after the fix (the fix
  only affects the previously-broken optimized-build codepath, not
  program logic).
- Release/RelWithDebInfo builds now run correctly on realistic synthetic
  datasets (see below) instead of crashing/misparsing.

This was necessary groundwork, not scope creep: without it there was no
way to build a working optimized `fastprot` to profile or later
benchmark against.

## Findings confirming/refining plan.md's "Target and scope"

- `Sequence` (`src/c++/Sequence.hpp`): confirmed, `std::string seq`, no
  protein-specific subclass, shared with DNA.
- `getAAInd()` (`ProtSeqUtils.cpp`): confirmed, maps the 20 canonical
  amino acids (case-insensitive via `toupper()`), returns sentinel `100`
  for anything else (gap, X/B/Z, whitespace, ...).
- `count_replacements()`: confirmed, calls `getAAInd()` on both
  characters and only tallies into the 20x20 `Matrix` when both are
  valid (`c1 != 100 && c2 != 100`) — unknown symbols silently excluded,
  as plan.md described.
- **Correction to plan.md**: `count_id_dist()` does **not** call
  `getAAInd()` at all. It does `toupper(*it1) == toupper(*it2)` directly
  on every position and divides by `s1.seq.size()`. This means:
  - Gaps and ambiguity codes are **not** excluded here the way they are
    in `count_replacements()` — a gap character participates in the
    comparison like any other symbol (two gaps at the same position
    count as a match; a gap vs. a residue counts as a mismatch).
  - There is **no length-equality check** between `s1` and `s2`; the loop
    walks `s1`'s iterator to `s1.seq.end()` while advancing `s2`'s
    iterator in lockstep. If `s2` is shorter, `it2` walks past
    `s2.seq.end()` — undefined behavior on real mismatched-length input,
    not a defined error path.
  - Phase 1/2's design must preserve *this* exact semantics (plain
    case-insensitive char compare, no alphabet filtering, denominator
    always `s1.seq.size()`) for `count_id_dist`'s replacement, which is
    different from `count_replacements()`'s semantics. Don't assume the
    two primitives share one "unknown symbol" behavior — they don't.
  - The undefined behavior on length-mismatched sequences should be
    treated as "existing behavior is UB / unspecified" — Phase 2 should
    pick a defined behavior (e.g. compare up to `min(len1,len2)`, or
    reject) and document it as a deliberate, disclosed change from "UB"
    to "defined", not as a silent behavior change to a case the old code
    ever legitimately handled.
- Build system / flags: confirmed. `ARCH_FLAGS = "-msse2"`
  (`src/c++/CMakeLists.txt:57`) is set but **dead** — its only consumer,
  `COMMON_FLAGS`, is commented out (line 204), so no `-msse2` (or any of
  the other flags in that commented block) is actually passed to the
  compiler today. There is no runtime SIMD dispatch anywhere in the
  codebase, confirming plan.md — any dispatch mechanism for the new
  protein code is new infrastructure.
- `RunExamples.sh` / `code_tests/`: confirmed as described — examples
  6/7 cover `fastprot`, diffed against `expected_output/`; `code_tests/`
  is not wired into the CMake build.

## Profiling results

Environment: Apple Silicon (arm64), AppleClang 15, `RelWithDebInfo`
build (`-O2 -g`, config.h fix applied), macOS `sample` (1ms sampling
interval). Synthetic dataset: 1500 protein sequences, 300 residues each,
5-40% divergence from a common ancestor (`prot_bench_big.fasta`,
generated for this audit — not checked in, regenerable).

**Kimura (`-D JCK`)**, 1500x300 dataset, wall time 1.23s for ~1.12M
pairs:
- 793/860 samples (92%) inside `calculate_kimura_dists` →
  `count_id_dist`.
- Within `count_id_dist` itself, the dominant cost is **`toupper()`**:
  roughly 690 of ~785 samples attributed to `count_id_dist` are actually
  inside `__toupper` calls (two per position: `toupper(*it1)` and
  `toupper(*it2)`), not the character comparison.
- The remaining ~66/860 samples are in output writing
  (`PhylipDmOutputStream::printPHYLIPfastSD` → `fwrite`), i.e. I/O, not
  a target for this work.

This is a more specific and more actionable finding than plan.md
anticipated: the Kimura path isn't just "dominated by the primitive" —
it's dominated by a **per-character libc function call**
(`toupper`) that a byte-encode-once-at-load-time design (already the
plan's Phase 1 approach) eliminates almost entirely, independent of and
in addition to whatever SIMD-compare speedup is layered on top. Encoding
residues to normalized (uppercase-collapsed) `uint8_t` codes at load
time removes the `toupper()` calls from the per-pair hot loop entirely,
since case normalization happens once per sequence instead of twice per
residue-pair-comparison.

**ED (`-D WAG`, non-ML) and ML (`-D WAG -m`)**: per user direction, not
profiled in detail this round — the ED path is not mature/not a current
priority. For context only: on the same 300x350 dataset used for initial
timing (before scaling up for Kimura profiling), wall time was ~2.0s for
ED and ~23.2s for ML, vs. 68ms for Kimura on identical input — i.e. ED is
~30x slower and ML ~340x slower than Kimura for the same pair count. This
is consistent with plan.md's expectation (and the "Future phases" section
on the Newton-Raphson/`expm()` solver) that ED/ML are dominated by the
numerical models consuming `count_replacements()`'s output, not by
`count_replacements()` itself — but this round did not profile ED/ML
directly to confirm the split quantitatively, so treat that as directional,
not verified.

## Bottom line for later phases

- Kimura path: profiling confirms the primitive (`count_id_dist`)
  dominates wall time, and specifically pinpoints `toupper()` overhead as
  the largest single cost within it — both the SIMD-compare work
  (Phase 2) and the byte-encode-at-load design (Phase 1) target this
  directly.
- `count_id_dist`'s actual semantics (no alphabet filtering, gaps
  participate normally, UB on length mismatch) differ from
  `count_replacements()`'s (filters via `getAAInd()`) — Phase 1's
  encode/decode design and Phase 4's differential tests need to account
  for both, separately, not assume one shared "ambiguity handling"
  story.
- ED/ML/`count_replacements()` deprioritized this round per direction
  received; revisit profiling them specifically if/when that path
  becomes a priority.
