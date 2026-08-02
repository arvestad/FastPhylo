# Phase 5 refactors: readability-function-cognitive-complexity

Before/after for each function refactored to address
`readability-function-cognitive-complexity` (28 findings, threshold 25).
Not every finding gets a full rewrite here - some are near-duplicate
functions where the same extraction is applied N times (documented
once, applied everywhere), and a few very close to the threshold are
left alone with a brief note on why splitting wouldn't actually help.

Every refactor below preserves behavior exactly - verified via full
rebuild + `ctest` + `RunExamples.sh` (byte-identical output) after each
one, same discipline as the rest of this lint pass.

## `fastdist/main.cpp`: `main()` (complexity 181 → all functions ≤ ~24)

By far the worst offender in the whole codebase - one 394-line `main()`
doing argument validation, translation-model setup, RNG seeding,
stream construction, and then three parallel output modes
(`--output-format=binary`, `--memory-efficient`, and the default
whole-matrix mode), each with its own no-bootstrap/bootstrap branching
and per-dataset loop.

**The key finding**: comparing the `--output-format=binary` block and
the `--memory-efficient` block side by side (via a normalized diff,
substituting the one literal `true`/`false` each used) showed they
were **structurally identical** - about 75 lines hand-duplicated,
differing only in a single boolean (`mem_eff_flag`) threaded through
`fillMatrixRow()`/`printRow()`. Nobody had parameterized it; they'd
just copy-pasted the block and changed `false` to `true`.

### Before (excerpt - the duplicated shape, ~150 lines total across both blocks)

```cpp
int main(int argc, char **argv){
  gengetopt_args_info args_info;
  TRY_EXCEPTION();
  sequence_translation_model trans_model;
  // ~80 lines: relaxng flags, model switch, tstvratio, bootstrapping
  // setup, srand(), ambiguity flags, fixfactor, ndatasets ...

  try {
    // ~35 lines: inputfilename/outputfilename resolution, istream/
    // ostream construction via two more switch statements ...

    if (args_info.output_format_arg == output_format_arg_binary) {
      StrFloRow dm;
      std::vector<Sequence> seqs;
      std::vector<string> names;
      std::vector<DNA_b128_String> b128seqs;
      Extrainfos extrainfos;
      for (int ds = 0; ds < ndatasets || ...; ds++) {
        std::string runId;
        if (!no_incl_orig && numboot == 0) {
          if (!istream->read(b128seqs, runId, names, extrainfos)) break;
          const size_t numberOfSequences = b128seqs.size();
          ostream->printStartRun(names, runId, extrainfos);
          ostream->printHeader(numberOfSequences);
          for (size_t i = 0; i < numberOfSequences; ++i) {
            fillMatrixRow(dm, b128seqs, trans_model, i, false);  // <-- literal
            dm.setIdentifier(names.at(i));
            if (useFixFactor) applyFixFactorRow(dm, fixfactor);
            ostream->printRow(dm, names.at(i), i, false);        // <-- literal
          }
        }
        else { /* ~35 more lines: bootstrap case, same shape */ }
        ostream->printEndRun();
      }
    } else if (args_info.memory_efficient_given) {
      // ~75 lines: THE SAME BLOCK AGAIN, byte-for-byte, except
      // every `false` above is `true` here.
    }
    else {
      // ~65 lines: a third, StrDblMatrix-based variant (fillMatrix()/
      // print() instead of fillMatrixRow()/printRow() - genuinely
      // different shape, not a duplicate of the other two).
    }
  }
  catch (...) { throw; }
  CATCH_EXCEPTION() catch (...) { cerr << "Unknown (non-Exception) error" << endl; }
  cmdline_parser_free(&args_info);
  return 0;
}
```

### After

`main()` itself is now argument parsing and dispatch only:

```cpp
int main(int argc, char **argv){
  gengetopt_args_info args_info;
  TRY_EXCEPTION();

  if (cmdline_parser(argc, argv, &args_info) != 0) exit(EXIT_FAILURE);

  handleEarlyExitFlags(args_info);

  sequence_translation_model trans_model = buildTranslationModel(args_info);

  int numboot = args_info.bootstraps_arg;
  bool no_incl_orig = args_info.no_incl_orig_given != 0u;
  seedRandom(args_info);

  bool useFixFactor = args_info.fixfactor_given != 0u;
  float fixfactor = args_info.fixfactor_arg;
  int ndatasets = args_info.number_of_runs_arg;

  try {
    std::unique_ptr<DataInputStream> istream = buildInputStream(args_info);
    std::unique_ptr<DataOutputStream> ostream = buildOutputStream(args_info);

    if (args_info.output_format_arg == output_format_arg_binary) {
      runRowStreamingOutput(*istream, *ostream, trans_model, ndatasets, numboot,
                             no_incl_orig, useFixFactor, fixfactor,
                             args_info.input_format_arg, RowOutputMode::Regular);
    } else if (args_info.memory_efficient_given != 0u) {
      runRowStreamingOutput(*istream, *ostream, trans_model, ndatasets, numboot,
                             no_incl_orig, useFixFactor, fixfactor,
                             args_info.input_format_arg, RowOutputMode::MemoryEfficient);
    } else {
      runFullMatrixOutput(*istream, *ostream, trans_model, ndatasets, numboot,
                           no_incl_orig, useFixFactor, fixfactor,
                           args_info.input_format_arg);
    }
  }
  catch (...) { throw; }
  CATCH_EXCEPTION() catch (...) { cerr << "Unknown (non-Exception) error" << endl; }
  cmdline_parser_free(&args_info);
  return 0;
}
```

`runRowStreamingOutput(..., RowOutputMode mode)` replaces both
duplicated blocks - called once with `RowOutputMode::Regular`, once
with `RowOutputMode::MemoryEfficient`. (First cut, called out in
review: this parameter was a bare `bool mem_eff_flag`, so the two call
sites read as `runRowStreamingOutput(..., false)` /
`runRowStreamingOutput(..., true)` - unreadable without checking the
signature. Replaced with the scoped enum shown above; the underlying
bool still gets reconstructed once, right at the point where the
pre-existing `fillMatrixRow()`/`printRow()` functions require it.) Its
own per-row inner loop was *itself* duplicated three times inside
(no-bootstrap case, bootstrap's "original sequences" case, bootstrap's
per-replicate case), so that became a further helper,
`fillAndPrintAllRows()`. Same pattern for the whole-matrix branch:
extracted to `runFullMatrixOutput()`, with its own three-times-duplicated
body factored into `setIdentifiersApplyFixFactorAndPrint()`. Six more
small helpers (`handleEarlyExitFlags`, `buildTranslationModel`,
`seedRandom`, `buildInputStream`, `buildOutputStream`) cover the
argument-handling preamble.

**A real bug risk found and avoided while doing this**: the first
extraction attempt moved `fillMatrix()` to *after* `printStartRun()` in
the no-bootstrap branch of the whole-matrix path - the original called
`fillMatrix()` first. Caught by re-diffing against the original source
before trusting the refactor, not by the test suite (both orders happen
to produce the same output here, since `printStartRun()` doesn't read
`dm`, so this wouldn't have shown up as a regression - but preserving
the exact original order removes the need to reason about that at all).

**Verification**: `RunExamples.sh` only exercises fastdist's plain
phylip/fasta/xml paths, none of `--output-format=binary`,
`--memory-efficient`, or `--bootstraps` - exactly the code this refactor
touched. Built the pre-refactor commit into a separate tree and diffed
outputs byte-for-byte against the refactored binary across six
combinations (binary, memory-efficient, bootstrapped phylip,
bootstrapped+binary, bootstrapped+memory-efficient, bootstrapped
+`--no-incl-orig`) - all identical, including stdout.

## `PhylipDmOutputStream::printPHYLIPfastSD` (complexity 86 → 0 findings in the file)

The batched phylip/XML matrix writer (shared by `fastdist`, `fastprot`,
and `XmlOutputStream`). All 86 points of complexity came from one
pattern repeated ~9 times: `if (writeXml || writeXmlSD) { ... } else {
... }`, re-checked at the header, at every row's prefix/suffix, and
three more times inside the innermost per-entry formatting loop
(not-finite case, large-value case, fast-path case).

**Before** (excerpt - the repeated-check shape, innermost loop only):

```cpp
void PhylipDmOutputStream::printPHYLIPfastSD(const StrDblMatrix &dm, FILE *out,
                                              bool writeXml, bool writeXmlSD,
                                              bool skipLeadingHeader) {
  ...
  for (size_t i = 0; i < numNodes; i++) {
    ...
    for (size_t j = 0; j < entriesPerRow; j++) {
      float f = dm.getDistance(i, j);
      if (!isfinite(f)) {
        if (writeXml || writeXmlSD) { row += "     <entry>-1</entry>\n"; }
        else { row += "        -1"; }
        continue;
      }
      ...
      if (intpart > 99) {
        if (f - (intpart * 1.0) < 0.000001) {
          if (writeXml || writeXmlSD) { n = snprintf(..., "     <entry>%10d</entry>\n", intpart); }
          else { n = snprintf(..., "%10d", intpart); }
        } else if (writeXml || writeXmlSD) { n = snprintf(..., "     <entry>%10f</entry>\n", f); }
        else { n = snprintf(..., "%10f", f); }
        ...
        continue;
      }
      ...
      if (writeXml || writeXmlSD) { row += "     <entry>"; row += &defstr[skip]; row += "</entry>\n"; }
      else { row.append(defstr.data(), 10); }
    }
    ...
  }
}
```

**After**: the two bools became one three-way enum (per the same
call-site-readability feedback as the `fastdist/main.cpp` refactor
above - `printPHYLIPfastSD(dm, fp, true, false)` at a call site forces
the reader to look up which bool means what; `printPHYLIPfastSD(dm, fp,
Format::Xml)` doesn't), and the entire per-entry formatting decision
tree - the actual source of the complexity - moved into three file-local
helpers (`appendEntry`, `appendFastEntry`, `appendRareEntry`) that
`printPHYLIPfastSD()` now just calls once per entry:

```cpp
enum class Format { Plain, Xml, XmlSD };  // in the class declaration

void PhylipDmOutputStream::printPHYLIPfastSD(const StrDblMatrix &dm, FILE *out,
                                              Format format, bool skipLeadingHeader) {
  const size_t numNodes = dm.getSize();
  const bool xml = format != Format::Plain;
  ...
  for (size_t i = 0; i < numNodes; i++) {
    row.clear();
    if (xml) { row += "    <row>\n"; entriesPerRow++; }
    else { /* name padding */ }

    for (size_t j = 0; j < entriesPerRow; j++) {
      appendEntry(row, dm.getDistance(i, j), format);
    }

    row += xml ? "    </row>\n" : "\n";
    fwrite(row.data(), sizeof(char), row.size(), out);
  }
  ...
}
```

`appendEntry()` holds the not-finite check and dispatches to
`appendFastEntry()` (the `|value| <= 99` fixed-width path,
`std::array<char, 11>`) or `appendRareEntry()` (the `snprintf()`
fallback for larger values) - each of those unconditional `if (writeXml
|| writeXmlSD)` branches now appears exactly once, in exactly one
helper, instead of being re-litigated at every call site.

Call sites updated to match: `print()`/`printSD()` in this same file
now pass `Format::Plain`/`Format::XmlSD`; `XmlOutputStream::print()`/
`printSD()` pass `Format::Xml`/`Format::XmlSD` (relying on
`skipLeadingHeader`'s existing default - XML output never sets it -
same as before this refactor).

**Out of scope, on purpose**: `fastprot_mpi/PhylipDmOutputStream.cpp`
has its own separate `printPHYLIPfastSD` with the same name and the
same bool pair - a different file, not touched here, consistent with
`fastprot_mpi` being out of scope for every phase of this lint pass.

**Verification**: full rebuild, `clang-tidy`'s cognitive-complexity
check back to zero findings for this file, `ctest`, and `RunExamples.sh`
(all 15 fixtures byte-identical) - this function is on the plain
phylip/xml output path that `RunExamples.sh` actually exercises
(unlike the fastdist refactor above), so no separate manual comparison
was needed this time.

## `computeDistance`/`computeTAMURANEIDistance` (complexity 65/85 → 0, plus `dist_level_4`'s 29 → 0)

The SIMD hot path: DNA_b128_String's b128-register-batched Jukes-Cantor
and Tamura-Nei distance functions (`computeDistance_DNA_b128_
String.cpp`, `computeTAMURANEIDistance_DNA_b128_String.cpp`). Handled
with more caution than any other Phase 5 target, for two reasons: it's
genuinely hot code (every pairwise sequence comparison in `fastdist`
goes through it), and - as this section explains - most of the
reported complexity here turned out not to be *our* code at all.

**First extraction, small payoff**: both functions open with a
fallthrough `switch` handling the ≤2 leftover b128s that don't divide
evenly into a `dist_level_1()` call (which always consumes 3 at a
time) - identical shape in both files, just with an extra `pyrts` term
in the Tamura-Nei version. Extracted verbatim into a
`sumRemainderBlocks()` helper in each file (kept as two separate
near-identical functions rather than merged, since the field counts
differ - the same call/duplicate-once tradeoff as the `PhylipDmOutput
Stream` section above). This only dropped complexity by 1 point each
(65→64, 85→84) - a `switch` counts as a single point in this cognitive-
complexity model regardless of how much code is inside its cases, so
moving a lot of code out doesn't move much complexity out with it. Worth
doing anyway (it's genuinely the least-readable part of either
function), but it wasn't the real story.

**The real finding**: re-running `clang-tidy` with full diagnostic
notes (not just the summary warning) showed almost the entire
complexity score - dozens of points - attributed to a single macro,
`shift_bytes_right_b128`, at every one of its call sites:

```
computeDistance_DNA_b128_String.cpp:364:3: note: +2, including nesting penalty of 1, nesting level increased to 2
  CONVERT_SUM_IMMEDIATEBYTESHIFT(total_sum_ts,total_sum_ts,EIGHT_BIT_MASK,1);
  ...
sse2_wrapper.h:233: note: expanded from macro 'shift_bytes_right_b128'
  #define shift_bytes_right_b128(__A, __IMM_BYTES)  _mm_srli_si128(__A,__IMM_BYTES)
simde/x86/sse2.h:1410: note: expanded from macro '_mm_srli_si128'
simde/x86/sse2.h:1381: note: expanded from macro 'simde_mm_bsrli_si128'
    if (HEDLEY_UNLIKELY(imm8 > 15)) { ... } else { ... }
```

On this machine (Apple Silicon), `_mm_srli_si128` isn't a native
intrinsic - it's simde's NEON-backed emulation (from `b5cef42`, "Fix
SSE2-on-ARM build via simde"), and simde's emulation macro contains a
runtime `if (imm8 > 15) {...} else {...}` to validate the shift amount
generically. Because `shift_bytes_right_b128` is a macro, that whole
if/else is *textually substituted* into `computeDistance()`/
`computeTAMURANEIDistance()` at every call site (there are over a
dozen between the two functions), and `clang-tidy`'s cognitive-
complexity check walks the expanded AST - so it scored every one of
those textually-inlined if/else pairs as if we had written that
branching by hand, nesting penalty and all. In practice `imm8` is
always a compile-time literal (`1`, `2`, `4`, or `8`) here, so the
compiler folds the check away entirely at `-O3` - it was never real
branching in the compiled output, only in what the linter saw.

**Fix**: converted `shift_bytes_right_b128`/`shift_bytes_left_b128`
from macros to templates in `sse2_wrapper.h` -
`shift_bytes_right_b128<N>(a)` instead of `shift_bytes_right_b128(a,
N)`. A non-type template parameter satisfies the same "must be a
compile-time immediate" requirement the macro existed for (documented
in both .cpp files' comments already, for the sibling `SUM_WITH_
PREVIOUS_LEVEL_IMMEDIATEBYTESHIFT`/`CONVERT_SUM_IMMEDIATEBYTESHIFT`
macros - a function parameter can't do this, `_mm_srli_si128` needs a
literal). The difference that matters here: a template function has
its own AST that a tool analyzes on its own terms, not text pasted
into every caller - so simde's internal if/else is now scored against
`shift_bytes_right_b128<N>()` itself (trivial), not against whichever
function happens to call it. Had to move the templates outside the
header's `extern "C" { ... }` block (this header is also compiled as
plain C - `sse2_wrapper.c.o` is a real build target - and templates
can't have C linkage), with a comment explaining why they're split out.

This dropped `computeDistance()` from 64 to 0 findings, `computeTAMURANEIDistance()`
from 84 to 0, and - as a side effect, since it uses the same macro -
`dist_level_4()`'s unrelated 29-point finding to 0 as well, without
touching `dist_level_4()` itself at all.

**Verification, given this is the hottest of hot paths**: `RunExamples.sh`'s
DNA fixtures (ex1/ex2/ex3) don't stress this enough to catch a subtle
regression, so: built the pre-refactor commit into a separate tree
(`build_old/`) and, on a synthetic 500-sequence/30,000bp dataset,
confirmed byte-identical `-O binary` output for both the default
(Jukes-Cantor) and `-D TN` (Tamura-Nei) distance functions, then timed
3 runs of each on old vs. new (`/usr/bin/time`) - within noise of each
other (~0.33s/~0.41s both before and after), consistent with the
change being AST-shape-only with no effect on the code the optimizer
actually emits. Also ran the full `ctest` + `RunExamples.sh` (15/15
byte-identical) suite as usual.

## `fnj`'s `XmlInputStream::readDM` (complexity 82 → 0)

`fnj`'s XML distance-matrix reader: a single 172-line function driving
a libxml2 pull-parser loop, structured as ten sequential `if (right
depth/parent-state/name) { switch(node_type) {...} }` blocks - one per
XML element it recognizes (`entry`, `row`, `dm`, `identity`,
`extrainfo`, `dms`, `run`, `identities`, `runs`, `root`). Each
condition chains 3-6 `&&`-ed booleans plus a depth check plus a name
comparison, and cognitive complexity charges for every one of those
operators, on top of the switch itself - ten times over.

**Key structural fact that made the extraction safe**: each block's
`name` comparison is mutually exclusive with every other block's (a
given XML node has exactly one tag name), so at most one block's full
condition can ever be true for a given node. That means "fall through
this `if` because the condition was false" and "handle it, then
`continue` the loop" are the only two things that ever happen, and
they're externally indistinguishable - both just mean "go read the
next XML node." That equivalence is what makes a clean split possible:
each block becomes a standalone method that re-checks its own guard
condition and returns `std::nullopt` if it doesn't apply, or optionally
a `readstatus` to return from `readDM()` immediately (`DM_READ` /
`END_OF_RUN` / `END_OF_RUNS` - the three cases that used to `return`
straight out of the middle of the function):

```cpp
std::optional<readstatus>
XmlInputStream::handleDmNode(int depth, int type, const xmlChar *name,
                              std::vector<std::string> &names, StrDblMatrix &dm) {
  if (!(l.in_root && l.in_runs && l.in_run && l.in_dms &&
        depth == 4 && xmlStrEqual(name, reinterpret_cast<const xmlChar *>("dm")) != 0)) {
    return std::nullopt;
  }
  switch (type) {
  case XML_READER_TYPE_ELEMENT:
    dm.resize(dmSize);
    l.in_dm = true;
    l.row_nr = -1;
    return std::nullopt;
  case XML_READER_TYPE_END_ELEMENT:
    l.in_dm = false;
    for (size_t namei = 0; namei < names.size(); namei++) { dm.setIdentifier(namei, names[namei]); }
    return DM_READ;
  default:
    return std::nullopt; // other node types (whitespace, comments, ...) intentionally ignored
  }
}
```

`readDM()` becomes pure dispatch - ten calls, each an `if (auto status
= handleXNode(...)) { return *status; }`:

```cpp
readstatus XmlInputStream::readDM(StrDblMatrix &dm, std::vector<std::string> &names,
                                   std::string &runId, Extrainfos &extrainfos) {
  int nr_of_ids = 0;
  int ret;
  while ((ret = xmlTextReaderRead(reader)) == 1) {
    if (xmlTextReaderIsValid(reader) != 1) { THROW_EXCEPTION("xml input does not validate"); exit(EXIT_FAILURE); }
    int depth = xmlTextReaderDepth(reader);
    int type = xmlTextReaderNodeType(reader);
    const xmlChar *name = xmlTextReaderConstName(reader);

    if (auto status = handleEntryNode(depth, type, name, dm)) { return *status; }
    if (auto status = handleRowNode(depth, type, name)) { return *status; }
    if (auto status = handleDmNode(depth, type, name, names, dm)) { return *status; }
    if (auto status = handleIdentityNode(depth, type, name, names, extrainfos, nr_of_ids)) { return *status; }
    if (auto status = handleExtrainfoNode(depth, type, name, extrainfos)) { return *status; }
    if (auto status = handleDmsNode(depth, type, name)) { return *status; }
    if (auto status = handleRunNode(depth, type, name, runId)) { return *status; }
    if (auto status = handleIdentitiesNode(depth, type, name, names, extrainfos)) { return *status; }
    if (auto status = handleRunsNode(depth, type, name)) { return *status; }
    if (auto status = handleRootNode(depth, type, name)) { return *status; }
  }
  return ERROR;
}
```

The ten handlers are private members (declared in `XmlInputStream.hpp`)
so they can read/write `l` (the `locator_t` parser-state struct) and
`dmSize` directly, same as the original inline code did.

Two genuinely dead locals (`const xmlChar *value;` and `bool run_read
= false;`, both declared but never referenced anywhere in the
function) were dropped while moving this code - confirmed unused via
grep before removing, not assumed.

**Verification**: full rebuild, `clang-tidy` back to zero findings,
`ctest`, `RunExamples.sh` (15/15 byte-identical - ex8/ex9 exercise
this function via `-I phylip`/`-I xml`). `RunExamples.sh` doesn't cover
every path through the ten handlers, so also built the pre-refactor
commit into a separate tree and manually diffed old vs. new on
`-I xml`, `-I phylip` with `-r 2 -O xml` (exercises the multi-run and
`identity`/`extrainfo` paths together), `--analyze-run-number`, and
the bad-input-count error path - all byte-identical.

## `fnj/main.cpp`: `main()` (complexity 69 → 0)

Refactored the same session as the `readDM()` split above, same file
family. `main()`'s complexity came from the usual mix (argument
validation, stream construction, a run loop) plus one genuine bug-shaped
finding: the run loop's body was branched on `args_info.input_format_arg
== input_format_arg_binary`, but **both branches of that if/else were
byte-for-byte identical** - `StrDblMatrix dm; for (...
istream->readDM(dm, ...) ...) {...}`, the exact same ~15 lines, copy-pasted
into both the `if` and the `else`. The input format only matters to
*which* `DataInputStream` implementation got constructed earlier
(`buildInputStream()`); by the time this loop runs, `istream->readDM()`
already dispatches virtually, so the branch was dead weight - probably
a leftover from when the binary path used a different matrix type,
later unified without ever removing the now-pointless split.

Extracted the shared body as `processRuns()` and deleted the branch
entirely (not "collapsed" - there was nothing left to collapse once
the duplicate was recognized as duplicate). Also extracted
`handleEarlyExitFlags()` (every early-exit check except the `isatty()`
one, which has to run before `cmdline_parser()` populates `args_info`
and so stays inline), `resolveMethods()`, `buildInputStream()`, and
`buildOutputStream()` - same shape as the `fastdist/main.cpp` helpers
from earlier in this document. Four more dead locals (`seqs`,
`b128seqs`, `latestReadSuccessful`, `speciesnames` - declared, never
read) were dropped, again confirmed unused via grep first.

**Verification**: full rebuild, `clang-tidy` back to zero findings,
`ctest`, `RunExamples.sh` (15/15 byte-identical - ex8/ex9 exercise
`main()` via `-r 2 -I phylip`/`-I xml`). Manually diffed old vs. new
on the same set of flag combinations as the `readDM()` verification
above (they share a binary and a test session), plus `--print-counts`,
`-m BIONJ`, `--print-relaxng-input`, and the too-many-input-files error
path - all byte-identical, including the now-branchless run loop
actually being exercised via both `-I binary` and `-I phylip`/`-I xml`
runs.

## `SequenceTree::mapSequencesOntoTree(istream&)` (complexity 61 → 0)

The three-overload family of `mapSequencesOntoTree` in `SequenceTree.cpp`
- `(char**, int)`, `(vector<Sequence>&)`, `(istream&)` - only the third
was flagged; the first two are ~15-line hash-map-and-lookup loops, the
third is a hand-rolled PHYLIP-ish sequence parser reading directly off
a `std::istream`, character by character.

Two duplicated "read one line of sequence characters" loops were the
main contributor: the first pass (each sequence's opening line, read
right after its name) and the "interleaved continuation" pass (later
lines of the same multi-line alignment, read once every sequence has
at least one line) both inlined the identical ~15-line
translate-and-append-until-newline loop. Extracted as a file-local
`readSequenceLine(istream&, string*)`. That alone only got the
function from 61 to 27 (still 2 over threshold) - the interleaved-
continuation pass's *surrounding* loop (`while (...) { for (...) {
peek/skip-name-field/skip-newlines; readSequenceLine(...); } }`) still
had real nesting of its own, so it became a second helper,
`readInterleavedContinuationLines()`, taking the same handful of
locals (`numSequences`, `seqlen`, `garbage`, `sequences`,
`actualNodeString`) the loop already closed over.

**No test coverage exists for this function** - it's called from
`Simulator.cpp` (simulation-tool support, not exercised by
`RunExamples.sh` or any wired-up test), and the only file that
exercises it directly (`tests/Tree_test.cpp`) isn't in any
`CMakeLists.txt`, so it never actually builds or runs. Wrote a small
throwaway standalone program instead (compiled against
`libfastphylo.a` directly, not committed) that builds a 3-leaf tree,
feeds a hand-crafted multi-line/interleaved input including an
unmatched name (exercises the "discard into `garbage`" path), and
dumps the result via `printSequencesPhylip()`. Built and ran this
against both the pre- and post-refactor library - byte-identical
output.

## `fastprot/main.cpp`: `main()` (complexity 59 → 0)

Same shape as the `fastdist`/`fnj` `main()` refactors earlier in this
document: `handleEarlyExitFlags()`, `buildTranslationModel()` (returns
`prot_sequence_translation_model` by value, same pattern as
`fastdist`'s `sequence_translation_model` - a small POD, safe to
return by value), `buildInputStream()`, `buildOutputStream()`, and a
`processRuns()` for the per-dataset/bootstrap loop.

**A real bug found and fixed**: the original order was

```cpp
gengetopt_args_info args_info;
TRY_EXCEPTION();
prot_sequence_translation_model trans_model;
#ifndef WITH_LIBXML
  if (args_info.input_format_arg == input_format_arg_xml){ ... exit ... }
#endif
if (cmdline_parser(argc, argv, &args_info) != 0) { exit(EXIT_FAILURE); }
```

- the `WITH_LIBXML=OFF` check reads `args_info.input_format_arg`
*before* `cmdline_parser()` populates it, i.e. reads an uninitialized
stack variable. In a `WITH_LIBXML=OFF` build, this is undefined
behavior that could spuriously exit with "built without XML support"
for a run that never asked for XML input, depending on whatever
garbage happened to be on the stack. Inert in the default build
(`WITH_LIBXML=ON` by `CMakeLists.txt`, so the `#ifndef` block doesn't
even compile in), which is presumably why it survived - `fastdist`'s
and `fnj`'s `main()`s both correctly call `cmdline_parser()` first.
Fixed by moving the check into `handleEarlyExitFlags()`, called after
`cmdline_parser()`, matching the other two.

`processRuns()` itself was still 3 over threshold after the first
pass (the original-data and bootstrap-replicate blocks each inlined
an identical "sd ? calculate-with-sd : calculate-without-sd; setIdentifiers"
step and an identical "print; sd-and-not-binary ? printSD : nothing"
step, but weren't otherwise identical - the original-data block also
runs `printStartRun()`/`printHeader()` in between, the bootstrap block
doesn't). Split those two common steps out as `calculateDistances()`
and `printDistances()`, called from both blocks with the
non-shared `printStartRun()`/`printHeader()` calls sandwiched between
them at the one call site that needs them.

**Verification**: full rebuild, `clang-tidy` back to zero findings,
`ctest`, `RunExamples.sh` (15/15 byte-identical). For the `WITH_LIBXML`
fix specifically: built a separate `WITH_LIBXML=OFF` tree and ran
`fastprot -I fasta` five times (undefined behavior isn't
deterministically reproducible, so this doesn't *prove* the old code
could fail, but confirms the new code never spuriously exits) plus
`fastprot -I xml` once, to confirm the "not built with libxml" error
still fires correctly when actually asked for.

## `fastphylo/io/XmlInputStream.cpp`: `XmlSequenceReader::readSequences` (complexity 55 → 0)

The shared XML sequence-file reader behind both `fastdist -I xml` and
`fastprot -I xml` (see this file's header comment on the Layout Phase C
merge). Same shape and same technique as `fnj/XmlInputStream.cpp`'s
`readDM()` split earlier in this document: five sequential "if (right
depth/parents/name) {...}" blocks (`seq`, `extrainfo`, `root`, `runs`,
`run`), each mutually exclusive by name, extracted into five
`handleXNode()` methods returning `std::optional<bool>` (nullopt = keep
reading; `true` = what the original returned inline from `run`'s
`XML_READER_TYPE_END_ELEMENT` case - the only place this function ever
returns early). `readSequences()` itself is now five dispatch calls.
Also dropped two dead locals (`value`, `run_read`) - same pair, same
non-use, as `fnj/XmlInputStream.cpp`'s `readDM()`, suggesting both were
written by copy-pasting the same template.

**Noted, not touched**: `XmlInputStream.hpp` already declares an unused
`enum streamstatus { RUN_NOT_FINISHED, RUN_FINISHED }` that looks like
an abandoned attempt to replace this function's `bool` return with
something more descriptive - exactly the enum-over-bool preference this
whole document has been applying. Didn't wire it up here: `bool
readSequences(...)` is a pure-virtual `DataInputStream` interface
implemented by several classes across both `fastdist` and `fastprot`
(`FastaInputStream`, `PhylipMaInputStream`, `XmlInputStream` in each),
all consumed via `if (!readSequences(...)) { break; }` - changing the
return type is a real, separate API-surface change touching a dozen-plus
files well beyond this one function's complexity fix, not something to
do as a drive-by.

**Verification**: full rebuild, `clang-tidy` back to zero findings,
`ctest`, `RunExamples.sh` (15/15 byte-identical - ex3's `fastdist -I
xml -O xml seq.xml` exercises this exact function).

## Task #47: remaining findings in the 26-47 range

The last batch of `readability-function-cognitive-complexity` findings
(~19, scattered across `include/fastphylo`/`src/fastphylo`/`src/c++/
apps/{fastdist,fnj,fastprot}`). Two of the twenty raw `clang-tidy`
hits (`tests/DNA_b128_String_test.cpp`, `tests/Tree_test.cpp`, both
`main()`) are skipped here: neither file is referenced by any
`CMakeLists.txt` (confirmed by grep), so neither builds - same
orphaned-test situation as `tests/Tree_test.cpp`'s other use noted
under the `mapSequencesOntoTree` section above. Refactoring code that
can't be compiled can't be verified, so these are left alone; whether
to wire them into the build at all is a separate question, not a
Phase 5 scope decision.

### `DNA_b128_String.cpp`: the four `correctDistanceWithAmbiguitiesUsing*()` overloads (26/26/26/26 → 0/0, 2 NOLINT'd at 26)

All four (`UsingBackgroundFrequences`/`UsingTransitionProbabilities`,
each with a `simple_string_distance` and a `TN_string_distance`
overload) turned out to run the *exact same* merge-style traversal
over `s1`'s and `s2`'s sorted ambiguity-position lists - only which
`compute_ambiguity_distance*()` function gets called, and which fields
of the result struct get updated, differs. Extracted the shared
traversal into two member function templates (one per
`compute_ambiguity_distance*()` family - `Impl` suffixed), parameterized
over the distance type and an `accumulate` lambda supplied by each of
the four thin public wrappers:

```cpp
template <typename Distance, typename Accumulate>
Distance DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequencesImpl(
    Distance sp, const DNA_b128_String &s1, const DNA_b128_String &s2, Accumulate accumulate) {
  Distance real_distance = sp;
  auto i1 = s1.ambiguities.begin();
  auto i2 = s2.ambiguities.begin();
  int pos1 = (i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
  int pos2 = (i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);
  while (pos1 < INT_MAX || pos2 < INT_MAX) {
    if (pos1 == pos2) {
      accumulate(real_distance, compute_ambiguity_distance((*i1).ambiguity, (*i2).ambiguity));
      ++i1; ++i2;
      pos1 = (i1 != s1.ambiguities.end() ? (*i1).position : INT_MAX);
      pos2 = (i2 != s2.ambiguities.end() ? (*i2).position : INT_MAX);
    } else if (pos1 < pos2) { /* ... */ }
    else { /* pos1 > pos2, symmetric */ }
  }
  return real_distance;
}

simple_string_distance
DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(simple_string_distance sp,
                                                                         const DNA_b128_String &s1,
                                                                         const DNA_b128_String &s2){
  return correctDistanceWithAmbiguitiesUsingBackgroundFrequencesImpl(sp, s1, s2,
      [](simple_string_distance &d, ambiguity_distance amdist) {
        d.transitions += amdist.purine_transition_prob + amdist.pyrimidine_transition_prob;
        d.transversions += amdist.transversion_prob;
        d.deletedPositions -= 1;
      });
}
```

Had to be a **member** function template, not a free function in an
anonymous namespace (the pattern used everywhere else in this
document) - it reads `s1.ambiguities`/`s2.ambiguities`, which is
private. Declared in `DNA_b128_String.hpp` (template member functions
need their declaration visible to callers even when, as here, the
definition stays in the `.cpp` and is only ever instantiated from that
same file).

This dropped all four *public* overloads' complexity to 0 (they're now
one-line forwarding calls), but the two `Impl` templates themselves
still measure 26 - one point over threshold. Left `NOLINT`'d rather
than split further: each is a single merge-style traversal over two
sorted lists with three mutually-exclusive cases
(`pos1==pos2`/`pos1<pos2`/`pos1>pos2`), already about as small as each
case can get; splitting the three cases into their own functions would
mean threading `pos1`/`pos2`/`i1`/`i2` through all of them by
reference, fragmenting one coherent loop for no readability gain -
exactly the "a few very close to the threshold are left alone with a
brief note" case `lint_plan.md` anticipated for this category.

**Verification**: full rebuild, `clang-tidy` down from 4 findings to
2 (both `NOLINT`'d), `ctest`, `RunExamples.sh` (15/15 byte-identical -
though none of those fixtures' sequences actually contain ambiguity
codes, so they only exercise the "empty ambiguity list, loop runs zero
times" path). Built a small phylip file with IUPAC ambiguity codes
(`N`, `R`, `M`) and compared the pre-refactor commit's `fastdist`
against the refactored one across `-D TN93`/`-D K2P` (the
`TN_string_distance`/`simple_string_distance` overloads
respectively) × with/without `-t` (`UsingTransitionProbabilities` vs
`UsingBackgroundFrequences`) × `-A`/`-R` (ambiguity-handling disabled,
and the pre-resolve step disabled) - all six combinations
byte-identical, covering all four refactored functions on data that
actually exercises their real logic, not just the empty-list path.

### `LeastSquaresFit.cpp`: `computeLeastSquaresEdgeLengths`/`computeLeastFloatSquaresEdgeLengths` (37/37 → 0/0)

Another double/float twin pair, both running the identical UNJ
edge-length-fitting algorithm - one on `StrDblMatrix`/`double`, one on
`StrFloMatrix`/`float`. Unlike the `DNA_b128_String.cpp` pair earlier,
the matrix types here don't share a common base with the same method
names for everything: `SequenceTree::tree2distanceMatrix()`/
`tree2FloatdistanceMatrix()` and the free functions `computeL2()`/
`computeFloatL2()` are differently *named*, not overloaded, so those
two call sites are passed into the shared template as function
pointers rather than resolved by overload resolution:

```cpp
template <typename Matrix, typename Real>
Real computeLeastSquaresEdgeLengthsImpl(const Matrix &orig_dm, SequenceTree &tree,
                                         void (SequenceTree::*treeToMatrix)(Matrix &),
                                         Real (*computeL2Fn)(const Matrix &, const Matrix &)){
  ...
  Matrix treeM(tree.getNumLeafs());
  (tree.*treeToMatrix)(treeM);
  return computeL2Fn(treeM, orig_dm);
}

double
computeLeastSquaresEdgeLengths(const StrDblMatrix &orig_dm, SequenceTree &tree){
  return computeLeastSquaresEdgeLengthsImpl<StrDblMatrix, double>(orig_dm, tree, &SequenceTree::tree2distanceMatrix, &computeL2);
}
```

The `Real` template parameter (not just `Matrix`) matters for
bit-for-bit fidelity: the original float version explicitly computed
`w1`/`w2`/`distChild1Child2` as `float` (narrowing from the `double`
arithmetic on the right-hand side at each assignment), not `double`
throughout - using `auto` instead of an explicit `Real` template
parameter would have silently changed every one of those roundings.

Consolidating the two 135-line bodies into one dropped both public
functions to 0, but left the shared `Impl` at 37 - still needed a
second pass, extracting four more helpers, each a self-contained,
previously-inline step of the per-merge bottom-up loop:
`buildRowIndexMaps()` (name↔row-index setup, was a leaf-check `if`
plus a not-found `if`+`USER_ERROR`), `computeWeightedNeighborSum()`
(the `sum +=` loop feeding both `EDGE(child1)`/`EDGE(child2)`),
`swapChildToLastRow()` (the row-bookkeeping `if` before
`removeLastRow()`), and `updateDistancesToParent()` (the final
per-row `setDistance()` loop). That got both `Impl` instantiations to
0.

**Also noted, not touched**: this pair's only callers are in
`src/c++/simulated_phylogenies/` (`BootstrapStats.cpp`/
`BootstrapTest.cpp`), which - like `tests/Tree_test.cpp` noted
earlier - aren't referenced by any `CMakeLists.txt`, so neither
`fastdist`, `fnj`, nor `fastprot` ever call this code. Out of scope to
fix (a build-wiring question, not a Phase 5 one), but relevant to how
this was verified.

**Verification**: full rebuild, `clang-tidy` back to zero findings,
`ctest`, `RunExamples.sh` (15/15 byte-identical, though - per the
above - this confirms no regression elsewhere, not this function,
which nothing in `RunExamples.sh` reaches). Wrote a standalone test
program (compiled against `libfastphylo.a` directly, not committed):
a 5-leaf tree (`((a,b),c,(d,e))`, chosen so the bottom-up loop
actually processes two non-root internal-node merges, exercising all
four new helpers rather than only the trivial 3-leaf-star case) with a
hand-built distance matrix, calling both `computeLeastSquaresEdgeLengths()`
and `computeLeastFloatSquaresEdgeLengths()` and dumping the L2 score
and every edge length at full precision. Built and ran against both
the pre- and post-refactor library - byte-identical output.

### `Sequences2DistanceMatrix.cpp`: `bootstrapSequences` (30 → 0), `DNA_b128_StringsFromPHYLIP` (42 → 0), `fillMatrix_Hamming`/`_JC`/`_K2P`/`_TN93` (31/31/41/47 → NOLINT'd)

`bootstrapSequences()` had an `if (32 < seqlen) {...} else {...}` split
that, on close reading, wasn't doing what it looked like. Both
branches called `rand()` to fill `samplePositions[0..seqlen-1]`; the
`32<seqlen` branch just did it in stride-32 batches (looks like an
unrolling attempt) while the `else` branch did it in one plain loop.
Since `rand()` is a serial call with no reordering, both visit
positions `0..seqlen-1` in the same order and call `rand()` the same
number of times - the batching changes nothing observable. Proved this
and unified both branches to the plain loop. The second half of each
branch (copying sampled positions into a bounded buffer and
`append()`-ing, chunked for sequences longer than the 16383-char
buffer) is *not* equivalent busywork - it's genuinely needed for
sequences that don't fit `buff` in one piece - but the `else` branch's
"just copy it all in one go" turns out to be the exact same code path
as the chunked version when `BUFFSIZE == seqlen` (which is always true
in the `else` branch's domain, since `stride=32 < 16383`): the chunk
loop runs zero iterations and the tail loop alone copies everything,
one `append()` call - identical to what the `else` branch wrote out by
hand. Unified onto the chunked version alone, deleting the `if/else`
and the entire `else` branch.

`DNA_b128_StringsFromPHYLIP()` had a real duplicate (the "has any
sequence reached `seqlen` chars" scan, computed once before the
interleaved-read loop and again inside it after every line) -
extracted as `anySequenceComplete()`. Preserved its literal semantics
("any", not "every" sequence complete) even though "every" reads more
intuitively for a stopping condition - not confident enough about the
intended PHYLIP-variant-tolerance behavior here to silently change it
mid-complexity-fix. The remaining interleaved-read `while` loop was
extracted as `readInterleavedLines()`, matching the technique used for
`SequenceTree::mapSequencesOntoTree(istream&)` earlier in this
document (also legacy hand-rolled multi-line PHYLIP parsing).

`fillMatrix_Hamming`/`_JC`/`_K2P`/`_TN93` were **not** touched - they're
four of the eight functions a separate, already-written project plan
(consolidating `Sequences2DistanceMatrix.cpp`'s duplicated
`fillMatrix*` family, not started) explicitly identifies as needing a
dedicated consolidation, ~90% identical to
each other and to their `fillMatrixRow_*` counterparts further down
this same file. A partial fix here would either duplicate or conflict
with that plan's eventual work, so left `NOLINT`'d with a comment
pointing at it instead.

**Verification**: full rebuild, `clang-tidy` down from 6 findings to 0
(4 `NOLINT`'d), `ctest`, `RunExamples.sh` (15/15 byte-identical - ex1/
ex2 exercise `DNA_b128_StringsFromPHYLIP()` directly). `bootstrapSequences()`'s
`rand()`-based resampling needed determinism to verify byte-for-byte,
so built the pre-refactor commit into a separate tree and compared
`fastdist -I phylip ... -b N -s <seed>` (fixed seed) output across
several bootstrap counts/seeds, `--no-incl-orig`, and - to actually
exercise the buffer-chunking path this file's `BUFFSIZE` logic exists
for - a synthetic 5-sequence/50,000-char dataset (`examples/seq.phylip`'s
sequences are only ~13-20 chars, short enough to never leave the
"else" branch that was just deleted) - all combinations byte-identical.

### `fnj/BinaryInputStream.cpp`: `readDM(StrFloMatrix&, ...)` (27 → 0)

Mechanical split: extracted the identifier-table-reading loop
(`readIdentifiers()`, a `while(true)` char-by-char scan for each
`:`-terminated name) and the distance-value-reading loop
(`readDistanceValues()`, the upper-triangular double loop with an
early-return-on-short-read) out of `readDM()`, no logic changes.

**This function turned out to be dead code** - `BinaryInputStream`'s
*other* overload, `readDM(StrDblMatrix&, ...)`, is the one `fnj`'s
`DataInputStream` interface actually declares `virtual`/`override`,
and that overload is still an unimplemented stub (a known, separate
gap - "needs writing + round-trip test" - not part of this lint pass's
scope). `fnj -I binary` always resolves to the stub and exits with
"Not implemented!" before this `StrFloMatrix` overload could ever run
- confirmed by grepping for callers (none) and by running `fnj -I
binary /dev/null` before and after this change (same stub message
either time).

**Verification**: given the above, there's no live code path to
exercise this function against - verification here is full rebuild +
`clang-tidy` (0 findings) + `ctest` + `RunExamples.sh` (15/15,
confirming no regression *elsewhere*) plus careful manual review that
the extraction is a pure mechanical move (every original statement
preserved verbatim, nothing reordered or altered) - not a claim that
the function's own behavior was runtime-verified, since it can't be
reached to test.

### `Sequence.cpp`: `readSequences` (47 → 0), `bootstrapSequences` (32 → 0)

A fourth (`readSequences`) copy of the same legacy interleaved-PHYLIP-
reading pattern already refactored twice earlier in this document
(`SequenceTree::mapSequencesOntoTree(istream&)`, `Sequences2Distance
Matrix.cpp`'s `DNA_b128_StringsFromPHYLIP()`) - same `anySequenceComplete()`
+ `readInterleavedLines()` split, same author ("Mehmood's Changes"
comment survives in all versions), same care to preserve this copy's
own small difference from the other two (`while (whileTrue)` here, vs.
`while (whileTrue || fin.eof())` in the other two - not homogenized,
since each copy may already encode its own subtly-different bug or
intentional variation that isn't mine to silently "fix" while chasing
a complexity number).

`bootstrapSequences` is the second copy of the `stride`-chunking
no-op-batching pattern this document already found and removed once
(`Sequences2DistanceMatrix.cpp`'s `bootstrapSequences()`).
`Sequence.cpp`'s own is actually simpler than that one: no
`BUFFSIZE`-bounded buffer here at all - it writes straight into a
pre-sized `std::string`, so unlike the `Sequences2DistanceMatrix.cpp`
version there wasn't even a genuine "long sequence" case to preserve.
The entire `if (stride < seqlen) {...} else {...}` collapsed to two
plain loops.

**Both turned out to be effectively dead code** in the shipped
binaries: `readSequences()` *is* reachable (`fastprot -I phylip`
routes through `PhylipMaInputStream`, which calls it - `RunExamples.sh`'s
`ex6`/`ex11`-`ex18` exercise it directly), but `bootstrapSequences()`'s
only real caller is `src/c++/simulated_phylogenies/BootstrapStats.cpp`
(same out-of-scope, unbuilt tool as earlier `LeastSquaresFit.cpp`
findings) - both `fastdist/main.cpp` and `buildtree.cpp` have it
present only as a commented-out call. `fastdist`'s actual bootstrap
path goes through the sibling `Sequences2DistanceMatrix.cpp` version
(`DNA_b128_String`-based, already verified earlier in this document).

**Verification**: full rebuild, `clang-tidy` back to zero findings in
both, `ctest`, `RunExamples.sh` (15/15 byte-identical, and - per
above - this really does cover `readSequences()`). For
`bootstrapSequences()`, since it can't be reached through any shipped
binary, wrote a small standalone program (compiled against
`libfastphylo.a` directly, not committed): three short seeded
sequences, five repeated bootstrap draws, output compared byte-for-byte
between the pre- and post-refactor library - identical.

### `SequenceTree_MostParsimonious.cpp`: `computeMostParsimoniousSequences` (45 → 0)

A Fitch-parsimony bottom-up DP: for each internal node and each
sequence position, computes the cost of assigning each of 4 symbols to
that position from its children's already-computed costs. That was a
triple-nested loop (symbol × child × child's-symbol) inside the outer
node/position loops - extracted as `computeParsimonyAtPosition()`,
leaving the caller with just the node/position loop nest:

```cpp
void computeParsimonyAtPosition(ParsimonyTree::Node *node, size_t j) {
  constexpr size_t numSymbols = 4;
  p_info &p = node->data[j];
  for ( size_t sym = 0 ; sym < numSymbols ; sym++ ){
    p[sym] = 0;
    ParsimonyTree::Node *child = node->getRightMostChild();
    for ( ; child != nullptr ; child = child->getLeftSibling() ){
      p_info &child_p = (child->data)[j];
      size_t score = child_p[0] + (sym != 0 ? 1 : 0);
      for ( size_t symC = 1 ; symC < numSymbols ; symC++ ){
        size_t symCp = child_p[symC] + (sym != symC ? 1 : 0);
        score = MIN(score, symCp);
      }
      p[sym] += score;
    }
  }
}
```

Left the leaf-initialization loop and the root/backtrack best-symbol
logic alone - one extraction was enough to clear the threshold, and
both of those touch `oldRoot->data.s.seq` as a shared mutable scratch
buffer in a way that looked riskier to pull out than the DP core,
which is a clean, self-contained computation with no such aliasing.

Like several other findings in this task, this function's only real
callers are in `src/c++/simulated_phylogenies/` (`BootstrapTest.cpp`/
`BootstrapStats.cpp`) - unreachable from `fastdist`/`fnj`/`fastprot`.

**Verification**: full rebuild, `clang-tidy` back to zero findings,
`ctest`, `RunExamples.sh` (15/15, confirming no regression elsewhere).
Wrote a standalone test (against `libfastphylo.a`, not committed): a
5-leaf tree (`((a,b),c,(d,e))`) with hand-picked sequences differing
at a few positions (enough to force real tie-breaks in the DP, not a
trivially-identical-sequences case), comparing the total parsimony
score and every node's reconstructed sequence between the pre- and
post-refactor library - byte-identical.
