# Plan: Linting and code-smell pass

## Goal

Bring the codebase to a clean, fully-linted state using the
`.clang-tidy` config already established in `modernization_plan.md`'s
Phase 0, without changing observable behavior except where a genuine
bug is found and disclosed. Same discipline as every prior phase of
this engagement: verify byte-identical output (`RunExamples.sh`,
`ctest`) after every step.

Branch: `lint-cleanup` (off `modernize-cpp17`'s tip - Lasse's call, a
new branch this time rather than continuing on `modernize-cpp17`
directly).

Scope: everything actually compiled into `fastphylo`/`fastdist`/`fnj`/
`fastprot` - i.e. `include/fastphylo/`, `src/fastphylo/`,
`src/c++/apps/{fastdist,fnj,fastprot}/`, `tests/`,
`benchmarks/bench_primitives.cpp`. `fastprot_mpi` and the smaller
untouched programs (`buildtree.cpp`, `CreateSimulatedData.cpp`,
`sequence_nj.cpp`) stay out of scope, consistent with every prior
phase's scoping. Confirmed-dead code (see Phase 0) is deleted, not
linted.

## Audit (done 2026-08-01)

Ran the existing `.clang-tidy` config across all 62 actually-compiled
files in scope. **11,257 findings total.** Top categories by volume:

| Category | Count | Nature (sampled, not assumed) |
|---|---|---|
| `bugprone-narrowing-conversions` | 363 | Mixed - needs its own triage, see Phase 2 |
| `readability-braces-around-statements` | 269 | Pure style |
| `readability-math-missing-parentheses` | 131 | Pure style |
| `readability-implicit-bool-conversion` | 101 | Pure style |
| `modernize-avoid-c-style-cast` | 58 | Same careful per-site treatment as modernization_plan.md's Phase 1, now needed in `dna/` (never got that pass) |
| `modernize-avoid-c-arrays` | 52 | Real type-signature changes |
| `readability-isolate-declaration` | 37 | Safe, mechanical |
| `modernize-loop-convert` | 34 | Index-based loops; same selective caution as modernization_plan.md Phase 3 |
| `readability-redundant-parentheses` | 31 | Safe, mechanical |
| `readability-function-cognitive-complexity` | 30 | Needs actual refactoring, not a lint fix |
| `modernize-use-auto` | 28 | Selective, per established convention (don't obscure important types) |
| `modernize-deprecated-headers` | 24 | Safe, mechanical (`<math.h>` -> `<cmath>` etc.) |
| `bugprone-throwing-static-initialization` | 23 | Sampled: false positives (SSE2 wrappers/`exp()` not `noexcept`-annotated, not actually throwing) |
| `readability-inconsistent-declaration-parameter-name` | 21 | Safe, needs per-site alignment |
| `readability-container-size-empty` | 21 | Safe, mechanical |
| `readability-inconsistent-ifelse-braces` | 19 | Safe, mechanical |
| `bugprone-easily-swappable-parameters` | 16 | API design smell, needs redesign not a fix |
| `readability-container-data-pointer` | 15 | Safe, mechanical |
| `bugprone-unchecked-string-to-number-conversion` | 15 | Sampled: real robustness gap - `atof`/`atoi` silently return 0 on malformed PHYLIP/XML numeric fields instead of reporting an error |
| `bugprone-switch-missing-default-case` | 13 | Sampled: not bugs - `fnj`'s XML parsing switches deliberately ignore non-element node types by design |
| `modernize-return-braced-init-list` | 12 | Safe, mechanical |
| `readability-simplify-boolean-expr` | 10 | Safe, mechanical |
| `bugprone-string-integer-assignment` | 8 | Sampled: false positives (`str += stream.get()` is deliberate one-char-at-a-time appending), but worth an explicit cast for clarity |
| `bugprone-macro-parentheses` | 8 | Sampled: real macros (`SUM_WITH_PREVIOUS_LEVEL` etc. in `dna/`) - none need text substitution (no token-pasting/stringification), so converting to `inline` functions is better than parenthesizing: real type-checking instead of just papering over precedence risk |
| `readability-redundant-string-init` | 7 | Safe, mechanical |
| `modernize-use-nullptr` | 7 | Leftover in `dna/` (never got modernization_plan.md's Phase 1) |
| `readability-else-after-return` | 5 | Safe, mechanical |
| `bugprone-exception-escape` | 5 | Sampled: real, narrow gap - all four `main()`s only catch `Exception`, not `std::exception`/`...`, so a non-`Exception` throw (e.g. `bad_alloc`) would `std::terminate` instead of printing a clean error |
| `bugprone-command-processor` | 5 | Sampled: `Simulator.cpp`'s `system()` calls (invoking SEQGEN) - needs confirming no untrusted input reaches them, not necessarily a fix |
| `readability-redundant-string-cstr` | 4 | Safe, mechanical |
| `modernize-use-emplace` | 4 | Mostly safe, spot-check copy/move semantics per site |
| `modernize-pass-by-value` | 4 | Calling-convention change, case-by-case |
| `bugprone-reserved-identifier` | 4 | Sampled: `_secant_search`/`_LEAST_SIGNIFCANT_BIT` - leading-underscore names, technically reserved in the global namespace. Safe rename. |
| `readability-use-std-min-max` | 3 | Safe, mechanical |
| `modernize-use-using` | 2 | Leftover in `dna/` |
| `modernize-use-bool-literals` | 2 | Safe, mechanical |
| `bugprone-random-generator-seed` | 2 | Sampled: `srand()` for bootstrap/simulation seeding - not a security context for this tool, no action needed |
| `bugprone-implicit-widening-of-multiplication-result` | 2 | Needs a per-site overflow-risk check |
| (several 1-2 count categories) | ~5 | Minor, folded into whichever phase touches that file |

Also known from earlier phases, explicitly flagged as candidates for
this pass: `aml/` (entirely unbuilt, confirmed by grep), `DNA_b128/
new_comp.cpp` (unreferenced by the build), `sequence_likelihood/
{DistanceMethodMatrix,LikelihoodMatrix}.hpp` (referenced only by `aml/`),
docbook's `GET_TARGET_PROPERTY(...LOCATION)` (deprecated CMake API,
`BUILD_DOCBOOK` defaults off), `STATIC`'s self-built gengetopt/libxml2/
zlib toolchain (FTP-hosted 2009-2011-era versions, unverified current
availability).

## Phases

Ordered by value-to-risk ratio, same principle as every prior phase in
this engagement: mechanical and safe first, judgment-requiring work
later, structural/design work last. Each phase gets its own
full-rebuild + `ctest` + `RunExamples.sh` verification, likely split
into multiple commits the way modernization_plan.md's Phase 1 and
Phase 6 were (one commit per pattern-type or per closely-related group,
not one giant commit).

### Phase 0 - Delete confirmed-dead code

`aml/` (6 files, entirely commented out of `CMakeLists.txt`),
`DNA_b128/new_comp.cpp`, `sequence_likelihood/DistanceMethodMatrix.hpp`
and `LikelihoodMatrix.hpp` (referenced only by `aml/`, confirmed by
grep before deciding, not assumed). Also the `DEBUG_MAIN`/`DEBUG_L1`-`L4`
macros in `computeDistance_DNA_b128_String.cpp` - all `if(false){INP;}`,
permanently compiled out, dead in the same sense. Zero behavior change
since none of this is reachable from any built target - cheapest,
safest phase, worth doing first to shrink the surface area for
everything else.

### Phase 1 - Correctness-relevant `bugprone-*` findings

Per-category, not blanket: `noexcept`-annotate the SSE2 wrapper
functions in `sse2_wrapper.h` (resolves `throwing-static-initialization`
correctly rather than restructuring init order for a non-issue); add
explicit `static_cast<char>` where `bugprone-string-integer-assignment`
flagged deliberate char-at-a-time string building; add `default: break;`
(with a comment - "other XML node types intentionally ignored") where
`switch-missing-default-case` flagged `fnj`'s XML parsing; convert `dna/`'s `SUM_WITH_PREVIOUS_LEVEL`-style macros to `inline`
functions (none need actual text substitution - real type-checking is
strictly better than parenthesizing);
rename the 4 reserved (`_`-prefixed) identifiers; decide and act on
`atof`/`atoi`'s silent-failure gap (likely: switch to `strtod`/`strtol`
with explicit error reporting on the PHYLIP/XML numeric fields, matching
this codebase's existing `THROW_EXCEPTION` error-handling convention);
decide on the `main()` exception-escape gap (likely: a final `catch
(...)` in each `main()`, printing something before exiting, rather than
a silent `std::terminate`); verify no untrusted input reaches
`Simulator.cpp`'s `system()` calls and document that finding either way;
check the 2 `implicit-widening-of-multiplication-result` sites for real
overflow risk given actual input sizes.

### Phase 2 - `bugprone-narrowing-conversions` triage (363)

**Status (2026-08-01): mostly deferred to a dedicated follow-up plan,
`distance_matrix_refactor_plan.md`.** Sampling (as planned below)
found the bulk (300+) trace back to one root cause: `DistanceMatrix`/
`FloatDistanceMatrix`/`DistanceRow`'s `int`-typed index accessors
(`getDistance`/`setDistance`/`getIdentifier`/`setIdentifier`) receiving
`size_t` at nearly every call site. Widening those parameters isn't a
safe mechanical fix, though - tracing it turned up a live `size_t
row = -1;` unsigned-wraparound sentinel trick in `Sequences2Distance
Matrix.cpp`'s `fillMatrixRow_*` family (and two more `field = -1`
assignments of unconfirmed status in `DistanceRow`/`FloatDistanceMatrix`'s
own stream constructors) that a careless widening pass could silently
break. That, plus the `fillMatrix*`/`fillMatrixRow_*` family being 8
near-identical, genuinely-too-long functions crying out for
consolidation (spotted independently by Lasse), made this real design
work, not a lint-pass line item - see `distance_matrix_refactor_plan.md`
for the full writeup, inventory, and phased execution plan.

What's fixed here directly (small, safe, confirmed no signed-arithmetic
risk): `nucleotide2ambiguity_nucleotide()`/`DNA_b128_String::
calcAmbiguityProbabilities()`'s `int basefreqA/C/G/T` params, widened
to `size_t` at the source (only ever added/multiplied/divided, never
subtracted).

What's left after the follow-up plan lands (not yet characterized):
`compute_Tamura_Nei`/`compute_Tamura_Nei_fixratio`'s `int strlen`/
`numAs`/etc - these must stay `int` (they subtract:
`strlen - sd.deletedPositions`) and need per-call-site explicit casts
instead; whatever remains once the DistanceMatrix-family sites are
gone will get a fresh sample-and-characterize pass rather than
assuming it's all the same shape.

### Phase 3 - Safe mechanical modernize-* cleanup, extending into `dna/`

`modernize-use-nullptr` (7), `modernize-deprecated-headers` (24),
`modernize-use-using` (2), `modernize-use-bool-literals` (2),
`readability-redundant-string-init` (7), `readability-container-size-empty`
(21), `readability-redundant-string-cstr` (4) - all the same safe,
mechanical patterns modernization_plan.md's Phase 1 already established
a process for, just extended to `dna/` (`DNA_b128`/`distance_methods`/
`sequence_likelihood`), which was deliberately excluded from that
plan's Phases 1-3 and never got this treatment. Also
`modernize-avoid-c-style-cast` (58) - same per-site
static_cast/reinterpret_cast/const_cast categorization discipline as
before (not a blind sweep), concentrated in `dna/`.

### Phase 4 - Broader readability/style mechanical fixes

`readability-braces-around-statements` (269), `readability-math-missing-parentheses`
(131), `readability-implicit-bool-conversion` (101),
`readability-redundant-parentheses` (31), `readability-isolate-declaration`
(37), `readability-inconsistent-declaration-parameter-name` (21),
`readability-inconsistent-ifelse-braces` (19), `readability-container-data-pointer`
(15), `modernize-return-braced-init-list` (12), `readability-simplify-boolean-expr`
(10), `modernize-use-emplace` (4, spot-check copy/move semantics per
site), `readability-use-std-min-max` (3), `readability-uppercase-literal-suffix`
(2), `readability-static-accessed-through-instance` (1),
`modernize-loop-convert` (34, same selective "only where the iterator
itself isn't needed" caution as before), `modernize-use-auto` (28,
same selective "don't obscure an important type" caution as before).
Biggest diff footprint of any phase by file count, but each individual
change is low-risk - verification cadence matters more here than
anywhere else.

### Phase 5 - Structural findings (real engineering, not mechanical lint)

`readability-function-cognitive-complexity` (30) - identify the worst
offenders, extract helper functions where that's a genuine readability
win, not a mechanical split; `modernize-avoid-c-arrays` (52) - convert
to `std::array`, which changes type signatures, so needs call-site
verification per instance; `bugprone-easily-swappable-parameters` (16)
- likely addressed opportunistically alongside whatever else touches
each flagged function, not as a dedicated sweep (reordering/wrapping
parameters in a struct is an API change with its own risk calculus per
call site); `modernize-pass-by-value` (4) - case-by-case, only where
the calling code already always passes a temporary or is fine with the
copy.

### Phase 6 - Tooling/build-system smells (done)

Two items flagged during the layout restructuring as candidates for
this pass, both low-risk since neither is exercised by default builds:

- docbook's `GET_TARGET_PROPERTY(...LOCATION)` (deprecated since CMake
  3.0). Turned out the result (`FASTDIST_EXE`) was never actually
  referenced anywhere in `docs/` - the "modern equivalent"
  (`$<TARGET_FILE:...>`) only works as a generator expression at the
  point of use, and there was no point of use to convert, so the
  correct fix was deleting the dead line rather than "modernizing"
  something with zero effect. `BUILD_DOCBOOK=ON` still configures
  cleanly after the removal (full docbook build remains untestable
  end-to-end here - `BUILD_DOCBOOK` defaults off, no `xmlto` in this
  environment - but the configure-time syntax is verified).
- `STATIC`'s self-built gengetopt/libxml2/zlib toolchain (FTP mirrors,
  2009-2011-era versions) - discussed directly rather than guessed at.
  Decision: remove `STATIC` entirely. Cross-platform (macOS/Windows/
  Linux) release binaries are wanted, via GitHub Actions - but that's
  new work needing its own plan (likely `vcpkg` for Windows/portable-
  macOS vendoring of *current* libxml2/zlib, not a revival of the old
  FTP-fetch mechanism), not a reason to keep 2009-2011-era
  FTP-mirror-dependent code around as a stopgap. gengetopt's own
  self-build fallback (independent of `STATIC` - triggers whenever no
  system `gengetopt` is on `PATH`, regardless) was left in place; it
  goes away naturally whenever gengetopt itself is migrated away from
  (separately tracked, not part of this pass - see "what's explicitly
  out of scope" below). See `github_actions_release_builds_plan.md`
  for the replacement.

### Phase 7 - Final verification and sign-off (done)

Full clean rebuild, `ctest`, `RunExamples.sh`, plus a benchmark
re-run (same discipline as modernization_plan.md's Phase 4) given
Phase 2's narrowing-conversion fixes and Phase 5's structural changes
both touch code on hot paths (`Matrix`, `ProtSeqCompare`, DNA distance
computation).

- Fresh `rm -rf build` + `cmake -DCMAKE_BUILD_TYPE=Release` +
  full parallel rebuild: clean, only pre-existing/unrelated warnings
  (macOS `sprintf` deprecation in the out-of-scope `Simulator.cpp`, an
  ignored x86-tuning flag on ARM).
- `ctest --output-on-failure`: 2/2 passed.
- `RunExamples.sh`: all 15 generated fixtures byte-identical against
  `expected_output/`. The 17-vs-15 file count in `expected_output/` is
  pre-existing (`ex4.out`/`ex10.out` aren't currently generated by any
  example script) - not a gap introduced by this pass.
- `bench_primitives`: healthy speedups across all three tested lengths
  (50/300/2000), 4x-57x depending on primitive/length, no regressions.
- Two additional ad hoc sanity checks beyond the plan's literal
  minimum, both exercising code Phase 2/5 touched directly: DNA
  distance computation (`fastdist`, all four models - JC/K2P/TN93/
  HAMMING - on a synthetic 400-seq/20000bp dataset) and protein
  maximum-likelihood distance (`fastprot -D WAG -m`, the
  `Matrix`/`MaximumLikelihood` hot path, on a synthetic 60-seq/500aa
  dataset). Both fast (well under a second) and consistent across
  repeated runs.

All 7 phases of this plan are now complete. Three follow-up efforts
were identified and spun out as their own plans along the way, none
started: `distance_matrix_refactor_plan.md` (Phase 2),
`phylip_reader_consolidation_plan.md` (Phase 5),
`github_actions_release_builds_plan.md` (Phase 6).

## What's explicitly out of scope

- `fastprot_mpi` (deferred per project_layout_plan.md, separate future
  request).
- `buildtree.cpp`/`CreateSimulatedData.cpp`/`sequence_nj.cpp` (never
  in scope for any phase of this engagement).
- Migrating off `gengetopt` (Lasse's earlier question, CLI11
  recommended - orthogonal to linting, not bundled here).
- A C++ namespace wrap (explicitly deferred in project_layout_plan.md's
  Decisions section).
