# Phase 5 benchmark results

Reproduce with `benchmarks/run_benchmarks.sh` (writes
`results_primitives.csv` and `results_end_to_end.csv`, checked in
alongside this file).

## Platform

Apple M1 Pro, macOS 26.4.1 (Darwin 25.4.0), AppleClang 15, arm64.
The SIMD path is SSE2 intrinsics translated to NEON at compile time via
simde (see phase1_design.md's dispatch decision) — **this is not a real
x86 SSE2 run**. No AVX2/AVX-512 dispatch exists (deliberately, per
phase1_design.md) so there's nothing further to test there. The scalar
remainder path (for lengths not a multiple of 16) is exercised by the
length-50 and length-300 cases below (48 and 288 are the largest
multiples of 16 below those lengths); it isn't a separate benchmarked
configuration, just confirmed correct by the same tests that pass at
those lengths (code_tests/ProtSeqCompare_test.cpp covers lengths from 0
to 3001, including many non-multiples of 16). **Not tested: a genuine
x86-64 machine, and a build/CPU without any of SSE2/NEON** (none was
available in this environment) — flagging this gap rather than
implying broader coverage than what was actually run.

## Primitive-level (`bench_primitives`, 3 warm-up + 21 timed reps, median ± IQR shown as [Q1, Q3])

| primitive | length | old median (µs) | new median (µs) | speedup |
|---|---|---|---|---|
| count_id_dist | 50 (short) | 2.00 [1.96, 2.00] | 0.38 [0.38, 0.42] | 5.3x |
| count_id_dist | 300 (typical) | 11.54 [11.50, 11.58] | 1.63 [1.58, 1.63] | 7.1x |
| count_id_dist | 2000 (long) | 74.13 [74.00, 74.50] | 12.00 [11.92, 12.04] | 6.2x |
| count_replacements | 50 (short) | 8.58 [8.54, 8.63] | 6.17 [6.17, 6.25] | 1.4x |
| count_replacements | 300 (typical) | 22.46 [22.42, 22.50] | 7.83 [7.79, 7.88] | 2.9x |
| count_replacements | 2000 (long) | 127.83 [127.75, 128.29] | 21.00 [20.96, 21.04] | 6.1x |

`count_id_dist` (the Kimura path's primitive, this plan's main target)
lands at **5-7x** across all three lengths. `count_replacements` (not
wired into fastprot this round — see phase0_audit.md/plan.md scope
notes; benchmarked here anyway since Phase 2 already implemented it) is
markedly *less* helped at short lengths (1.4x at 50 residues) than long
ones (6.1x at 2000) — expected, since it's the scalar tally per
plan.md's design decision 3, not a SIMD path; its speedup here comes
entirely from encode-once-per-sequence removing `getAAInd()`'s
`toupper()` calls from the hot loop, not from vectorization, and that
fixed per-call overhead matters more at short lengths. Reported
honestly per plan.md's explicit ask not to obscure this with an
aggregate number.

## End-to-end (`fastprot -D JCK`, wall clock, 7 reps/case, median; `/dev/null` output to isolate compute from I/O)

| n_seqs | length | old (s) | new (s) | speedup |
|---|---|---|---|---|
| 150 | 100 | 0.067 | 0.038 | 1.8x |
| 150 | 500 | 0.206 | 0.060 | 3.4x |
| 600 | 300 | 1.657 | 0.304 | 5.5x |
| 1500 | 300 | 10.234 | 1.709 | 6.0x |

End-to-end speedup rises with problem size (1.8x at the smallest case,
approaching the primitive-level ~6-7x by 600-1500 sequences) because
fixed per-process overhead (startup, parsing, one-time model setup)
doesn't shrink with the comparison speedup, so it's proportionally more
of the total at small sizes. This matches expectations, not a red flag
— real workloads (dozens-to-thousands of sequences, per plan.md's Phase
4 discussion of realistic dataset sizes) sit well inside the
size range where the full primitive-level speedup shows through.

Note: these end-to-end numbers write to `/dev/null` specifically to
isolate the comparison-primitive speedup this plan targets from the
separate, already-logged output-writing bottleneck
(`PhylipDmOutputStream::printPHYLIPfastSD`'s fwrite()-per-entry cost —
see phase0_audit.md's "Post-Phase-3 finding"). Real invocations that
write a large distance matrix to disk will see a smaller end-to-end
win than shown here, dominated increasingly by that separate, not-yet-
fixed bottleneck as matrix size grows — same reasoning as the ED/ML
Newton-Raphson item, logged as follow-up work rather than folded in.

## Acceptance bar

plan.md asks for a target speedup factor decided with Lasse before this
phase; none was set before starting. Reporting the measured numbers
above for Lasse to judge against, rather than asserting a pass/fail
this round.
