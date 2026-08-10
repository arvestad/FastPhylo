# Phase 1 design — protein sequence byte encoding

Covers: the alphabet/code table, encode/decode functions, and the SIMD
dispatch strategy, per plan.md's Phase 1. Scope stays what Phase 0
confirmed: `fastprot`'s `count_id_dist()` (Kimura's primitive, no
alphabet filtering) and `count_replacements()` (20x20 tally, filters via
`getAAInd()`).

## Input alphabet, confirmed from the actual parsers

`FastaInputStream.cpp:74` rejects (throws) any character outside:

    abcdefghiklmnopqrstuvwyzxABCDEFGHIKLMNOPQRSTUVWYZX -.?

That's the 26 letters minus `j`, case-insensitive (25 letters), plus
space, hyphen, dot, question mark — 29 distinct symbol-identities after
case folding. FASTA input is guaranteed to stay within this set.

`PhylipMaInputStream.cpp` has **no such check** — it delegates to
`Sequence::readSequences`, which is the shared DNA/protein reader and
accepts arbitrary bytes. So a PHYLIP-format protein file can contain any
character; the current code's behavior on such input is whatever
`toupper()`+direct comparison happens to do (defined per-character, but
outside the "intended" alphabet) for `count_id_dist`, and silent
exclusion (via `getAAInd()`'s sentinel) for `count_replacements` if the
byte isn't one of the 20 canonical letters.

This means the encoding can't just cover "20 AAs + a couple of
ambiguity codes" and lump everything else into one bucket without
changing behavior for legal-but-unvalidated PHYLIP input. See "OTHER
bucket" below for how this is handled and what, precisely, changes.

## Code table

One 256-entry lookup table, `uint8_t residue_code[256]`, built once
(static init). Single source of truth, `ProtSeqCode.hpp`. Layout:

| Codes | Meaning |
|---|---|
| 0–19  | Canonical amino acids, **same order as `getAAInd()`** (`ProtSeqUtils.cpp`): A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V |
| 20–24 | The 5 other FASTA-legal letters that aren't canonical AAs: B,O,U,X,Z (26 letters minus `j` minus the 20 canonical = these 5) |
| 25    | `-` (gap/indel) |
| 26    | ` ` (space) |
| 27    | `.` (dot) |
| 28    | `?` (question mark) |
| 29    | OTHER — catch-all for any byte not covered above |

N = 30 codes, fits comfortably in a `uint8_t` (plan.md estimated "~25";
the real FASTA-legal alphabet is 29 symbols, so 30 codes including the
OTHER bucket). Both upper- and lower-case letters map to the same code
(`residue_code[c] == residue_code[toupper(c)] == residue_code[tolower(c)]`
for letters), matching `count_id_dist`'s and `getAAInd()`'s existing
case-insensitivity.

**One table serves both primitives:**
- Mismatch counting (replaces `count_id_dist`): two positions match iff
  their codes are equal. Because the table is built from
  `toupper()`-folded identity, this is **exactly** `count_id_dist`'s
  existing per-position comparison for every character FASTA input can
  ever contain (codes 0–28 are all injective over distinct
  characters/case-classes). No behavior change for FASTA input.
- Replacement tally (replaces `count_replacements`): valid iff
  `code < 20`. This is exactly `getAAInd()`'s canonical/non-canonical
  split restated as a range check — codes 20–29 are exactly the set
  `getAAInd()` maps to sentinel `100` today.

**OTHER bucket (code 29) — the one disclosed behavior change:** any byte
outside the 29-symbol FASTA-legal set collapses to a single code. Two
*different* out-of-alphabet bytes (e.g. a literal `#` and a literal `!`
from a hand-edited PHYLIP file) would now compare as a *match*, whereas
today's `toupper()`-based comparison would call them a mismatch (they're
different characters). This only affects PHYLIP input containing bytes
FASTA input could never legally contain in the first place — it's a
narrow, disclosed corner case, not a change to any behavior the
`RunExamples.sh` fixtures or realistic protein data would ever exercise.
Flagging this explicitly per design-decision 4 ("no behavior change to
public results... within a justified tolerance") — this is the
justification, and it should be called out again in Phase 4's test plan
and Phase 6 documentation, not silently shipped.

Sequence **length mismatch**: Phase 0 found `count_id_dist` has actual
UB here (walks `s2`'s iterator past `s2.end()` if shorter than `s1`).
The new primitive will compare up to `min(len1, len2)` and treat any
length difference as contributing to the mismatch count for the
remaining positions of the longer sequence, with the denominator staying
`len(s1)` (matching the old denominator choice) — turning UB into a
defined, sensible behavior rather than leaving it unspecified. This is a
new decision, not a preserved one; call it out again in Phase 4.

## Encode/decode functions (`ProtSeqCode.hpp`/`.cpp`)

- `uint8_t encode_residue(char c)` — table lookup, O(1).
- `void encode_sequence(const std::string &seq, std::vector<uint8_t> &out)`
  — encodes once per sequence at load time, not per comparison. This is
  what eliminates the `toupper()`-per-comparison cost Phase 0's profile
  identified as the dominant cost in `count_id_dist` (92% of Kimura's
  wall time; ~88% of that inside `count_id_dist` itself was `toupper`
  calls) — encoding folds case exactly once per residue per sequence,
  not twice per residue per *pair*.
- `char decode_residue(uint8_t code)` — for error messages / any future
  serialization of the new representation; not performance-sensitive.
  Canonical uppercase representative for codes 0–28; `'?'` for the OTHER
  code 29 (matches the existing FASTA-legal `?` = "unknown" convention).

## SIMD dispatch strategy

plan.md flags there's no existing runtime CPU-feature dispatch anywhere
in FastPhylo to reuse (Phase 0 confirmed this, and also found
`ARCH_FLAGS=-msse2` is dead code — not actually passed to the compiler
today).

Decision: **no runtime dispatch for this pass.** Byte-wise SIMD compare
only needs the SIMD baseline of each target ISA:
- x86-64: SSE2's `pcmpeqb`/`pmovmskb` is part of the x86-64 baseline —
  every x86-64 CPU has it, no runtime check needed.
- AArch64 (incl. Apple Silicon): NEON is part of the AArch64 baseline,
  same reasoning. `simde` (already a build dependency as of commit
  `b5cef42`, used for `DNA_b128`) translates the same `_mm_*` SSE2 call
  sites onto NEON at compile time — reuse that, so one intrinsics-based
  implementation, written once against SSE2 intrinsics, works on both
  ISAs without conditional code paths.
- Anything else (32-bit x86, other architectures without simde
  coverage): plain scalar loop, selected at compile time via the same
  `#if`/CMake platform check `DNA_b128`'s `sse2_wrapper.h` already uses.

This is a deliberate scope reduction from plan.md's "build multiple code
paths (scalar, SSE4, AVX2) with a runtime CPU-feature check" — AVX2
would need runtime dispatch (not all x86-64 CPUs have it), which is real
extra infrastructure for a further speedup on top of an already-large
win from (a) removing `toupper()` from the hot loop and (b) SSE2/NEON
baseline vectorization. Flagging this for Lasse: if AVX2 turns out to
matter after Phase 5's benchmark, adding runtime dispatch on top of this
design is additive, not a rewrite — the scalar and SSE2/NEON paths this
phase builds stay as the fallback and baseline tiers.

## What's reimplemented vs. untouched this pass

- New: `ProtSeqCode.hpp`/`.cpp` (alphabet table, encode/decode) — Phase 1,
  this commit.
- New (Phase 2): SIMD mismatch-count primitive operating on
  `std::vector<uint8_t>`, replacement-matrix tally primitive (scalar,
  per plan.md decision 3).
- Untouched this pass: `Sequence` (`std::string seq` stays the storage
  type; encoding happens as a derived/cached representation, not a
  replacement of `Sequence` itself — avoids touching the many other
  `Sequence` call sites (DNA, I/O, bootstrapping) outside this project's
  scope), `getAAInd()`, `count_id_dist()`, `count_replacements()`
  (old implementations stay in place per design-decision 5, reachable
  behind a switch added in Phase 3).
