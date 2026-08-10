# Plan: fix the binary distance-matrix format's run-boundary gap

## Goal

`fastdist`/`fastprot -O binary` and `fnj -I binary` (implemented in
`project_fnj_binary_input_gap`) round-trip correctly for the common
case - one dataset, any number of bootstrap replicates sharing one
header. But the on-disk format has **no marker for where one run's
data ends and the next run's header begins**. `fastdist -r N -O binary`
(multiple datasets) writes N full headers back-to-back into one
stream with nothing distinguishing "more matrix floats" from "a new
header" - confirmed by direct experiment
(`fastdist -I phylip <2-dataset file> -r 2 -O binary | fnj -I binary
-r 2`): `fnj`'s first run silently swallows dataset 2's header bytes
as if they were more distance floats (misreporting `count=3` where it
should be `1`), and the second run gets nothing. Exit code 0 - no
error, just wrong output.

This plan fixes it by giving the format explicit run framing, at the
cost of a breaking wire-format change. That's acceptable here: nothing
could actually read the format at all until two days ago (the reader
was a stub), so there are no real persisted files or external readers
to stay compatible with - the format is this project's own internal
detail.

## Why not fix the reader alone

Runs `original data + N bootstrap replicates` all share one header
(`printStartRun()`/`printHeader()` write once; `print()`/`printRow()`
are called once per replicate) - this is a deliberate compactness
choice already baked into the format (avoids rewriting the names table
for every replicate). Because of that, nothing in the byte stream
itself distinguishes "one more matrix body" from "a new header" -
both are just bytes, and a matrix body's raw floats can't be
told apart from a tag string without already knowing which one to
expect. The reader cannot recover a boundary that was never written.
The only real fix is for the **writer** to record, once per run, how
many matrix bodies follow - then the reader can count instead of
guess.

## The fix

**Wire format** (`BinaryDmOutputStream::printHeader()` /
`BinaryInputStream::readHeaderIfNeeded()`):
- Bump the tag from `"FASTPHYLO 1"` to `"FASTPHYLO 2"` (same 11-byte
  length, so `tagLength` stays unchanged) - signals the shape change
  and lets the reader reject a foreign/corrupt file with a clear error
  instead of silently misparsing it. This is also a real, disclosed,
  pre-existing gap worth closing alongside: today the tag is read into
  a variable and **never validated** at all.
- Add a `long matricesPerRun` field, written right after the node
  count and before the names table: how many `print()`/`printRow()`
  matrix bodies follow this header before the next header (or EOF).

**Writer** (`BinaryDmOutputStream`): add a `matricesPerRun` constructor
parameter (no change to the `DataOutputStream` virtual interface -
`printHeader()`'s signature and every other writer stay untouched).
Checked both real call sites: in both `fastprot` and `fastdist`,
`matricesPerRun = (opts.no_incl_orig ? 0 : 1) + opts.bootstraps` is a
single value fixed for the whole program invocation (not something
that could vary per dataset/run today), and `opts` is already in scope
at both `buildOutputStream()` definitions - so this is a same-function
change, no signature threading needed through `main()`.

**Reader** (`BinaryInputStream`): replace the single `bool
input_was_read` (which today can only ever fire once, for the entire
stream) with a small run-boundary state machine:

- `NeedHeader` - parse the next header (validate tag, read node count
  + `matricesPerRun` + names fresh). If the stream has nothing left at
  all, return `END_OF_RUN` (true end of file - same "no more real
  data" signal `PhylipDmInputStream` already uses, so `fnj`'s existing
  `-r`-bounded outer loop keeps working unchanged).
- `InRun` - reuse the cached names, read one more body, return
  `DM_READ`, decrement the remaining-bodies counter.
- `JustFinishedRun` - the call right after a run's last body was
  delivered: return `END_OF_RUN` without touching the stream yet (this
  is required by `processRuns()`'s `for (; readDM(...)==DM_READ; )`
  loop shape, which needs a dedicated call to learn "stop" before the
  next real read starts), then go back to `NeedHeader` for the call
  after that. This exactly mirrors how `XmlInputStream` already
  signals run boundaries via `END_OF_RUN`/`END_OF_RUNS` - `fnj`'s
  reading loop doesn't need to change at all.

Settled member-variable shape (replaces the current single `bool
input_was_read`):
```cpp
enum class RunState { NeedHeader, InRun, JustFinishedRun };
RunState state = RunState::NeedHeader;
size_t bodiesRemainingInRun = 0;
```
A malformed/truncated body read (short read where `bodiesRemainingInRun`
was still `> 1`, i.e. the file claimed more bodies than it actually
has) resets to `NeedHeader` and returns `END_OF_RUN` rather than
throwing - consistent with every other reader's existing "a short read
just means done" convention, not a new failure mode to invent.

## One dead-code item folded in, not deferred

**`readDM(StrFloMatrix&, ...)`** (the sibling overload) has zero
callers anywhere in the codebase - re-confirmed by a fresh repo-wide
grep during Phase A (the only other hits are its own declaration/
definition and `XmlInputStream`'s identically-dead, identically-named
stub, which is out of this plan's scope). Duplicating the new state
machine into genuinely dead code is pure waste. Deleting it outright
in this plan rather than keeping it in sync for its own sake - matches
this session's simulator-code-deletion precedent
(`project_simulator_code_removal`).

## fastprot_mpi: investigated, explicitly NOT touched here (revised from the original draft)

The original draft of this plan assumed `apps/fastprot_mpi/
BinaryDmOutputStream.{hpp,cpp}` (+ `PhylipDmOutputStream` sibling) were
simple orphaned leftovers - delete them, fix the one `#include`, done.
Phase A's direct-inspection requirement (not just grep) found that's
wrong, and the real situation is more entangled:

- `apps/fastprot_mpi/DataOutputStream.{hpp,cpp}` **is** genuinely
  compiled (it's in `FASTPROT_MPI_SRCS`) and has **diverged** from the
  shared `include/fastphylo/io/DataOutputStream.hpp` - different
  virtual signatures (e.g. `printRow(StrFloRow&, name, int row)` vs.
  the shared class's `printRow(StrFloRow&, name, size_t row, bool
  mem_eff_flag = false)`), no `printBootstrapSpliter`, pre-`#pragma
  once` style. `fastprot_mpi/main.cpp` declares its `DataOutputStream
  *ostream` against *this* local class.
- `apps/fastprot_mpi/BinaryDmOutputStream.hpp`/`PhylipDmOutputStream.hpp`
  **are** genuinely orphaned (not in `FASTPROT_MPI_SRCS`, their `.cpp`
  bodies are never compiled) - but they `#include "DataOutputStream.hpp"`
  unqualified, so they subclass fastprot_mpi's *local*, diverged
  `DataOutputStream`, not the shared one.
- Simply pointing their `#include`s at the shared `io/` headers (the
  original draft's plan) would make them subclass the *shared*
  `DataOutputStream` instead - a completely different, ABI-incompatible
  type from the `DataOutputStream *ostream` pointer `main.cpp` actually
  holds. That would not compile (assigning an unrelated type to
  `ostream`), not just fail to link.
- Properly reconciling this needs fastprot_mpi's own `DataOutputStream`/
  `XmlOutputStream` (its real, compiled, diverged base classes) migrated
  onto the shared hierarchy - which is exactly the "fold `fastprot_mpi`
  into `fastprot`" work `project_layout_plan` already scoped out as its
  own separate future request, not a side effect of a wire-format fix.

**Conclusion**: `fastprot_mpi`'s `-O binary`/`-O phylip` output paths
are *already* latently broken (an object of the wrong dynamic type
family assigned through a mismatched base-pointer, if it ever linked
at all) independently of this plan - `BUILD_WITH_MPI` defaults off and
is unbuildable in this environment, so nothing observes this today.
This plan's writer-side constructor change doesn't make a working
thing newly broken; it was never going to correctly build in the first
place. Left untouched here, disclosed instead - real fix belongs to
the deferred `fastprot_mpi` fold-in, not this plan. (Recorded as its
own memory note so it isn't lost before that fold-in happens.)

## Phases

- **A** - DONE (2026-08-07). Settled the design with no code changes:
  confirmed the wire format (tag bump, field placement), confirmed the
  reader state machine, re-confirmed `readDM(StrFloMatrix&, ...)`'s
  zero-caller status, re-confirmed `matricesPerRun = (no_incl_orig?0:1)
  + bootstraps` is never reassigned after computation in either
  `fastprot`'s or `fastdist`'s `main()` (grepped every assignment site).
  Direct inspection of the `fastprot_mpi` files (not just grep, per
  this phase's own requirement) overturned the original draft's
  "delete + fix the include" plan for it - see the revised section
  above. That single correction is exactly why this phase exists
  before any code gets written.
- **B** - DONE. Writer: `BinaryDmOutputStream` constructor +
  `printHeader()` changes; updated `fastprot`'s and `fastdist`'s
  `buildOutputStream()`. (`fastprot_mpi` intentionally not touched -
  see above.) Compiled clean.
- **C** - DONE. Reader: implemented the state machine in
  `BinaryInputStream::readDM(StrDblMatrix&, ...)`, added tag
  validation (throws on a mismatched/foreign tag instead of silently
  misparsing - a real, disclosed pre-existing gap, the tag used to be
  read and never checked at all); deleted the dead `StrFloMatrix`
  overload. Compiled clean. **Verified live, before touching any
  tests**, against the exact scenario found broken when this plan was
  written: `fastdist -I phylip <2-dataset file> -r 2 -O binary | fnj
  -I binary -r 2` now reports `count=1` for both runs (previously
  `count=3` for run 1 and an empty run 2). Also re-verified the two
  already-working cases (single dataset, bootstrap replicates sharing
  one header) still work unchanged.
- **D** - DONE. Extended `tests/BinaryDmIO_test.cpp` with
  `test_multiple_runs_are_correctly_separated()` (two full runs, each
  own header+body, different names/sizes per run). Updated the two
  existing cases for the new constructor parameter. **Caught a real
  test-design bug of my own while writing it**: first draft gave run 1
  and run 2 different `matricesPerRun` (1 vs. 2) on the same output
  stream instance - but that value is fixed per-instance by design
  (matching the real callers, where it never varies per dataset), so
  the writer silently wrote 1 unclaimed extra float body past what the
  header promised, and the assertion caught it immediately
  (`END_OF_RUN` where a second `DM_READ` was expected). Fixed by giving
  both runs the same `matricesPerRun` (1), which is what any real
  caller would do - not a reader bug, a wrong test.
- **E** - DONE. Added example 20:
  `fastdist -I phylip seq.phylip -r 2 -O binary | fnj -I binary -r 2
  -O xml -d 4` - discovered `seq.phylip` (example 1's fixture) already
  contains two concatenated datasets, so no new fixture file was
  needed. Both runs correctly report `count=1` (previously: `count=3`
  for run 1, empty run 2). Also had to regenerate `expected_output/
  ex16.out` (the direct byte-diff of `fastprot -O binary`'s raw
  output) - expected, since the wire format itself changed (tag +
  matricesPerRun field); confirmed its new bytes are exactly
  `"FASTPHYLO 2"` + node count + `matricesPerRun` before saving it as
  the new expected output, not just copied blindly. Full clean
  rebuild (0 errors) + `ctest` 4/4 + `RunExamples.sh` 20/20
  byte-identical + `RunCliChecks.sh` all green.
- **F** - DONE. This section itself is the sign-off; memory updated
  (see `project_fnj_binary_input_gap` and this file's own memory
  entry).

## Decisions

- Breaking wire-format change accepted outright, no migration path -
  justified above (no real consumers existed before this week).
- Both dead-code cleanups above (the `StrFloMatrix` overload,
  `fastprot_mpi`'s orphaned local files) are in scope for this plan
  rather than separate follow-ups, since leaving either would either
  waste effort keeping dead code in sync or introduce a new latent
  build trap as a direct side effect of this plan's own writer-side
  change - not pre-existing issues being opportunistically swept in.

## Status

Written 2026-08-07. All six phases done the same day. Verified: clean
rebuild, `ctest` 4/4, `RunExamples.sh` 20/20 byte-identical,
`RunCliChecks.sh` all green. Not yet committed.
