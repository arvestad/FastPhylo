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
