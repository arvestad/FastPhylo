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

