# Plan: migrate off gengetopt

## Goal

Replace gengetopt (the `.ggo` DSL, its generated C parser, and the
FTP-fetch-and-self-build fallback that provisions it when no system
`gengetopt` is on `PATH`) with a modern, header-only C++ CLI library.
Two explicit constraints from direct discussion:

1. **Preserve the current CLI interfaces as closely as possible** -
   same flags (long and short forms), same types, same defaults, same
   enum choice strings, same positional-`FILE` behavior, same
   early-exit flags (`--print-relaxng-*`). Existing scripts and
   muscle memory shouldn't break.
2. **Make adding/changing options easier going forward** - Lasse is a
   fan of Python's `argparse`; the replacement should let a new option
   be added with one fluent call at the point of use, not a separate
   DSL file plus a build-time code-generation step.

This is a sibling effort to [[project_lint_cleanup]] (which
deliberately left gengetopt alone - "Migrating off `gengetopt`...
orthogonal to linting, not bundled here") and directly de-risks one of
[[project_github_actions_release_builds]]'s open questions: that plan
flags Windows as needing `gengetopt` provisioned via `vcpkg` (no
system package exists there); a header-only replacement needs no
provisioning on any platform, so this migration should land *before*
that plan's Phase D (Windows).

## Current state, verified directly (not assumed)

Four apps use gengetopt, each via its own `apps/<name>/gengetopt/cmdline.ggo`:

| App | Options | Notes |
|---|---|---|
| `fastdist` | 18 | most enums, most `_given` checks (tstvratio/pyrtvratio/fixfactor/seed all optional-with-fallback) |
| `fastprot` | 14 | shares most of its option set with `fastprot_mpi` |
| `fnj` | 10 | smallest surface, good first-migration candidate |
| `fastprot_mpi` | 13 | near-duplicate of `fastprot`'s `.ggo` (already flagged as a duplicate elsewhere - see `project_layout_plan.md`); **out of scope for every phase so far**, same as the rest of this engagement |

Confirmed via grep: every reference to `args_info`/`cmdline_parser`/
`gengetopt_args_info` is confined to each app's own `main.cpp` - no
other `.cpp`/`.hpp` in the tree touches the generated struct directly.
That's a direct result of the lint-cleanup pass's cognitive-complexity
refactors (`lint_phase5_refactors.md`), which already extracted each
app's option-consuming logic into small helper functions:
`handleEarlyExitFlags`, `buildTranslationModel`/`buildInputStream`/
`buildOutputStream`/`seedRandom` (fastdist), `resolveMethods`/
`buildInputStream`/`buildOutputStream` (fnj), equivalents in
`fastprot`. Those helpers already take `args_info`-shaped state by
parameter and contain no parsing logic themselves - only the
option-parsing block and the specific `.foo_arg`/`.foo_given` field
reads need to change; control flow does not. This makes the migration
meaningfully lower-risk than it would have been before that refactor.

Confirmed live (`fastprot --help`/`--version`/a bad enum value), so the
"preserve interfaces" goal has a concrete baseline to check against:

```
$ fastprot --help
Usage: fastprot [OPTION]... [FILE]...
  -h, --help                    Print help and exit
      --detailed-help           Print help, including all details...
  -V, --version                 Print version and exit
  -o, --outfile=filename        output filename...
  -D, --distance-function=ENUM  Distance function  (possible values="ID",
                                  "JC", ... default=`WAG')
  ...

$ fastprot --version
fastprot 1.0.10

$ fastprot -D BOGUS
fastprot: invalid argument, "BOGUS", for option `--distance-function' (`-D')
$ echo $?
1
```

Every enumerated option (`input-format`, `output-format`,
`distance-function`, `ambiguity-frequency-model`, `method`, `speed`)
compiles to gengetopt's own generated C enum (e.g.
`distance_function_arg_JC`), switched on throughout each app's
`main.cpp`. The positional `FILE` argument is optional, single-valued
in every `.ggo` (`args_info.inputs_num`/`args_info.inputs[0]`) - no app
currently accepts more than one positional file, confirmed by reading
all four `.ggo` files in full, not assumed.

`RunExamples.sh` exercises real invocations (`fastdist -I phylip
seq.phylip`, `fnj -r 2 -I phylip dm.phylip`, `fastprot ... -D JCK -O
binary`, etc.) - short options, `-X value` spacing, enum values, and
flags are all covered by the regression suite. **`--help`, `--version`,
and error paths (bad enum, unknown flag, `--print-relaxng-*`) are
not** - those need manual verification during migration, not just a
green `RunExamples.sh`.

## Chosen library: CLI11

[CLI11](https://github.com/CLIUtils/CLI11) - header-only C++11,
BSD-3-Clause, explicitly designed around a fluent
`app.add_option(...)`/`app.add_flag(...)` call-per-option API that is
the closest C++ analog to `argparse.add_argument()` available today.
Recommended over the alternatives:

- **cxxopts**: also header-only, but modeled more on
  `boost::program_options`'s builder-chain style than argparse's
  imperative one - a real option, but less of a match for the stated
  preference.
- **Boost.Program_options**: adds a full Boost dependency this project
  doesn't otherwise have, for a capability a single header already
  covers - unjustified weight.
- **Hand-rolled parser**: rejected outright - defeats "make adding
  options easier," and this project already has a hand-rolled-parser
  cautionary tale in `arg_utils_ext.cpp`'s own ad hoc flag handling
  found during the lint pass.

CLI11 covers every feature the current `.ggo` files use: typed
options (`string`/`int`/`float`), enumerated choices via
`CLI::CheckedTransformer` or `CLI::IsMember` (validates *and* converts
in one step - a direct fit for converting a `--distance-function`
string straight into a project-defined `enum class`), flags, default
values (auto-shown in `--help`), an optional single positional
argument, and both `-o value`/`-o=value`/`--long=value` syntaxes.
`.count("opt") > 0` is the direct equivalent of gengetopt's
`..._given` fields used throughout (e.g. `seed_given`,
`fixfactor_given`, `analyze_run_number_given`).

**Vendoring**: recommend vendoring a single pinned CLI11 header
in-tree (`third_party/CLI11.hpp` or similar), matching this project's
existing pattern for `simde` - avoids adding a network dependency
(FetchContent) or a package-manager dependency (vcpkg/apt/brew) to a
clean build, and sidesteps exactly the kind of provisioning gap that
motivated this migration in the first place. Specific version to pin
is a Phase A decision (pick a recent stable release at implementation
time, not decided here).

## What this removes

Once all in-scope apps are migrated: every `.ggo` file, the
generated-code `add_custom_command`s and their `MAKE_DIRECTORY` calls,
`find_program(GENGETOPT gengetopt)` and its FTP-fetch-and-build
fallback (`build_gengetopt.sh`, the `ftp://ftp.gnu.org/...
gengetopt-2.22.2.tar.gz` download, `GENGETOPT_VERSION`) - all deleted
outright, not left as a NOLINT'd fallback "just in case" (same
convention as `STATIC`'s removal - see
[[feedback_modernize_dont_keep_legacy_stopgaps]]). This is the last
piece of 2009-2011-era vendored-toolchain machinery in the build.

## Decisions / constraints to hold to

- **Flag names, short letters, types, defaults, and enum choice
  strings are preserved exactly** for every option in scope. This is
  a hard constraint, not a "close enough."
- **Positional `FILE` argument stays optional and single-valued**,
  falling back to stdin when absent - matches every app's current
  behavior, verified above.
- **`--print-relaxng-input`/`--print-relaxng-output`'s early-exit
  behavior is preserved** - these must fire (and conflict-check
  against each other) before any other required-option validation,
  same as today's `handleEarlyExitFlags` helpers already isolate.
- **`--help`/`--version` text is explicitly allowed to look
  different.** CLI11's autogenerated help formatting differs
  structurally from gengetopt's; byte-matching it isn't the goal or a
  reasonable ask. What must be preserved is *content*: every option,
  its short/long name, and a description at least as clear as today's.
  Show a before/after side-by-side for the first migrated app (`fnj`)
  and get explicit sign-off before repeating the pattern on the rest,
  rather than assuming the same judgment call three more times.
- **Exit codes and stderr behavior on bad input should match**
  (gengetopt exits `1` with a one-line stderr message on an invalid
  enum value or unknown flag, confirmed above) - CLI11 does the same
  via its default `CLI11_PARSE` error handling; verify per app rather
  than assume.
- **No new CLI functionality added as part of this migration** -
  subcommands, config-file loading, shell completion, etc. are all
  things CLI11 *can* do and gengetopt couldn't, but adding them now
  would conflate "migrate" with "extend." The ease-of-adding-options
  win is about *future* work, not licence to expand scope here.

## Open questions - settled 2026-08-04

- **`fastprot_mpi` scope**: **in.** Leaving one straggler app on
  gengetopt would mean the old FTP-fallback machinery could never be
  fully deleted - direct instruction to include it.
- **Vendoring mechanism**: **vendored header**, confirmed - a single
  pinned `CLI11.hpp` in `third_party/`, not `FetchContent`/vcpkg/a
  system package. (Note: this project's one existing header-only
  dependency, `simde`, actually uses `find_path()` against a
  system/brew install, *not* in-tree vendoring - an earlier draft of
  this plan claimed otherwise. CLI11 is this project's first
  genuinely vendored-in-tree dependency, not a repeat of an existing
  pattern.)
- **Automated CLI-level tests**: **add a small script.**
  `examples/RunCliChecks.sh` (added in Phase B) checks, per migrated
  app: `--version` exit code and exact output, `--help` exit code and
  presence of every long option name, and that an invalid enum value
  and an unrecognized flag both exit nonzero. Wired into
  `.github/workflows/build-and-test.yml` right after
  `RunExamples.sh`. Extend it with one more `check_app` call per app
  as Phases C/D land.

## Phases

### Phase A - Settle the open questions above (done)

See "Open questions - settled 2026-08-04" above.

### Phase B - Vendor CLI11, prove the pattern on `fnj` (done)

Vendored `third_party/CLI11.hpp` (v2.7.2, the latest release at the
time, confirmed via the GitHub releases API - not guessed). Wired into
`src/c++/CMakeLists.txt` as a plain include path on `fnj`'s target
only (`target_include_directories(fnj PRIVATE ... third_party)`) -
`fastdist`/`fastprot`/`fastprot_mpi` are untouched and still build via
gengetopt, confirmed by a full build succeeding with both toolchains
present simultaneously.

Migrated `fnj/main.cpp`: `gengetopt_args_info` replaced by a plain
`FnjOptions` struct (values plus `_given` bools for the options whose
*presence* matters, matching gengetopt's `..._given` fields exactly);
`input-format`/`output-format` got their own scoped enums,
`--method` binds directly to `NeighborJoining.hpp`'s existing
`NJ_method` domain enum via `CLI::CheckedTransformer` (eliminating the
old `resolveMethods()` translation function entirely - CLI11 doing the
validation *is* the translation). `--dm-per-run` is accepted but still
unused, matching pre-migration behavior exactly (confirmed via grep
before touching it - nothing in `fnj` ever read
`args_info.dm_per_run_arg` either). Deleted `apps/fnj/gengetopt/`
(the `.ggo` and its directory) and every fnj-specific line in
`src/c++/CMakeLists.txt`'s generated-code plumbing.

Two real gaps found and fixed by diffing against the captured
baseline, not assumed away:
- **`--help` text leaked CLI11's internal validator repr** -
  `CheckedTransformer`'s default behavior appends `value in {a->0,
  b->1,...} OR {0,1,...}` to every enum option's description, which is
  strictly worse than gengetopt's clean `(possible values=...)`.
  Fixed with `.description("")` on each `CheckedTransformer`, since
  the possible-values/default text is already written by hand in each
  option's own description.
- **Exit codes didn't match.** gengetopt exits `1` uniformly on any
  parse problem (confirmed live: bad enum value, unrecognized flag -
  both `1`). CLI11 assigns distinct codes per error subtype (105 for
  a failed validator, 109 for an unexpected argument, etc.) and only
  `0` for `--help`/`--version`. Fixed by catching `CLI::ParseError`,
  calling `app.exit(e)` for its printed message, then explicitly
  `std::exit`-ing `0` if the error is a `CLI::Success` (help/version)
  or `1` otherwise - this is exactly the kind of gap the plan's own
  "verify per app, don't assume" note about error paths was for.

Verified: full clean rebuild (all four apps, two toolchains
coexisting) + `ctest` (2/2) + `RunExamples.sh` (byte-identical,
including fnj-exercising examples 8/9) + the new
`RunCliChecks.sh` (all checks pass) + a manual side-by-side of
`fnj --help`/`--version`/a bad enum value/an unknown flag against the
baseline captured pre-migration. **`--help` formatting sign-off is
pending** per this plan's own "get sign-off on the first migrated app"
rule - not yet repeated on Phase C/D until that lands.

### Phase C - Migrate `fastdist`

Most complex app: 18 options, the `sequence_translation_model`-building
logic (`tstvratio`/`pyrtvratio`/`no-tstvratio`/`fixfactor`, several
`_given` checks feeding conditional defaults), two enums
(`input-format`×3, `output-format`×3, `distance-function`×4,
`ambiguity-frequency-model`×2). Same verification discipline as Phase
B: `RunExamples.sh` (examples 1-3), manual help/version/error check,
plus the DNA-distance-model sanity checks already used in
`lint_plan.md`'s Phase 7 (all four `-D` values) since this is
hot-path-adjacent code.

### Phase D - Migrate `fastprot` (and `fastprot_mpi`, if Phase A includes it)

14 options, the largest single enum (`distance-function`×10). If
`fastprot_mpi` is in scope, do both together given how much of their
option set is identical - matches the "near-duplicate, migrate as a
pair" treatment already used for other `fastprot`/`fastprot_mpi`
duplication elsewhere in this project. `RunExamples.sh` covers
`fastprot` extensively (examples 6/7/11-18, including the ML path and
binary-output path); `fastprot_mpi` has no example coverage today
(confirmed absent from `RunExamples.sh`) - needs the same manual
verification discipline as `--help`/`--version` above, or a decision
to add minimal coverage first.

### Phase E - Remove all now-dead gengetopt infrastructure

Once every in-scope app is migrated: delete the remaining `.ggo`
files, `find_program(GENGETOPT)`, the FTP-fetch `add_custom_command`s,
`build_gengetopt.sh`, `GENGETOPT_VERSION`, and every `MAKE_DIRECTORY`/
generated-file plumbing in `src/c++/CMakeLists.txt` tied to them - full
removal, matching the `STATIC` precedent, not a kept-around fallback.

### Phase F - Final verification and sign-off

Full clean rebuild + `ctest` + `RunExamples.sh` (byte-identical, same
discipline as every prior phase in this engagement), plus the
manual `--help`/`--version`/error-path smoke-test matrix across every
migrated app (since none of that is covered by the automated
regression suite). Update `github_actions_release_builds_plan.md` to
note its Windows-gengetopt open question is now moot, once this lands
before that plan's Phase D.

## What's explicitly out of scope

- Any new CLI functionality (subcommands, config files, completion)
  beyond what gengetopt already provided.
- `fastprot_mpi`, unless Phase A explicitly decides otherwise.
- Re-litigating option names/behavior "while we're in there" - this is
  a mechanical migration, not a UX redesign. Any option Lasse actually
  wants changed should be its own follow-up, after this lands.
