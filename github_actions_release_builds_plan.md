# Plan: GitHub Actions cross-platform release builds

## Status

**A-D done (2026-08-07)**: `build-and-test.yml` now runs on Linux,
macOS, and Windows on every push/PR, all three green. Ten real,
disclosed portability/behavior bugs found and fixed along the way -
see Phase D below for the full list. **E-F not started**: artifact
packaging, the tag-triggered release job, and the version-source
mechanism (settled in Phase A, not yet implemented) remain.

## Goal

Automated GitHub Actions builds producing distributable binaries for
**macOS, Windows, and Linux**, for "reasonably modern OS versions" of
each, triggered on tagged releases (and ideally exercised on every
push/PR too, to catch breakage before a release rather than at it).

This builds on the existing `.github/workflows/build-and-test.yml`
(Linux-only build+test, no release artifacts - see "What already
works" below), not a from-scratch design.

## There's real prior art, and a likely diagnosis of why it stalled

A `.github/workflows/release.yml` already existed on this branch's
history and was deleted:

- `8a905a0` added `build-and-test.yml` (the Linux-only job that
  survives today).
- Separately, `release.yml` (a 3-OS matrix release workflow) was
  present earlier, got a fixup commit **`ec8ae49` "Adding complexities
  for OpenBLAS"**, then was deleted entirely in **`6ac573b` "Do not
  have time for problem-solving github actions"**.

Retrieved from git history (`git show 6ac573b^:.github/workflows/release.yml`),
the deleted workflow:

- Matrixed `ubuntu-latest`/`macos-latest`/`windows-latest`.
- Installed BLAS per OS: `apt-get install libopenblas-dev
  liblapack-dev` (Linux), `brew install openblas` (macOS, forcing
  `BLA_VENDOR=OpenBLAS` and pointing `CMAKE_PREFIX_PATH` at Homebrew's
  prefix), `microsoft/vcpkg-action` installing `openblas` (Windows,
  with `CMAKE_TOOLCHAIN_FILE` pointing at vcpkg).
- **Never installed libxml2 or gengetopt, on any OS.**
- Packaged binaries into a zip/tar.gz per OS and uploaded them via
  `softprops/action-gh-release` on `v*` tags.

The commit sequence (add release.yml → "adding complexities for
OpenBLAS" → delete it, giving up) points at two compounding problems:

1. **BLAS/LAPACK cross-platform detection was already painful enough
   to need a fixup commit before the workflow was abandoned.** Forcing
   `BLA_VENDOR=OpenBLAS` uniformly on every OS (rather than each
   platform's natural choice) is extra complexity for no clear payoff
   here - the project's actual `CMakeLists.txt` just does plain
   `find_package(BLAS REQUIRED)`/`find_package(LAPACK REQUIRED)` with
   no vendor hint, and that already resolves correctly to Apple's
   Accelerate framework on macOS with zero extra configuration
   (confirmed repeatedly this session, every local build). There's no
   evidence in the codebase that OpenBLAS specifically is required
   anywhere.
2. **libxml2 (and gengetopt) were never provisioned on Windows at
   all**, where neither exists by default - unlike Linux/macOS, which
   at least have *something* find-able (apt package, system dylib) even
   if imperfect. `WITH_LIBXML` defaults `ON`, so `find_package(LibXml2)`
   would very likely have failed outright on the Windows runner. This
   looks like the more fundamental gap - the BLAS complications may
   have been a red herring the previous attempt chased while the real
   blocker (no libxml2 on Windows) went undiagnosed.

Worth registering as the working hypothesis for Phase A below, not a
certainty - nobody has re-run that workflow to confirm exactly where it
failed.

## What already works (the foundation, not a redo)

`.github/workflows/build-and-test.yml`, currently on `master`/
`modernize-cpp17`/`speed2026a` push and every PR:

```yaml
runs-on: ubuntu-latest
# apt-get install cmake libxml2-dev liblapack-dev libblas-dev
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
# RunExamples.sh, diff every fixture against expected_output/
```

This proves two things worth relying on rather than re-litigating:

- The **current** default CMake configuration (post `STATIC` removal
  in `lint_plan.md`'s Phase 6 - no forced BLAS vendor, plain
  `find_package(LibXml2)`) genuinely works end-to-end on a real
  GitHub-hosted Ubuntu runner, not just on this dev's local macOS
  machine.
- `RunExamples.sh` is already wired up as the project's regression
  smoke test and is directly reusable as each new platform's
  build-verification step.

Note: this job does **not** `apt-get install gengetopt` - meaning
every run today exercises the FTP-based gengetopt self-build fallback
(`find_program(GENGETOPT gengetopt)` failing, then downloading
`gengetopt-2.22.2.tar.gz` from `ftp.gnu.org` and building it). Kept
deliberately when `STATIC` was removed, since gengetopt migration is
tracked separately - but it means this CI job currently depends on an
FTP mirror's availability on every single run, worth fixing
independently of the rest of this plan (see Phase B).

## The three platforms' real dependency gaps

- **Linux** (`ubuntu-latest`): solved by the existing job. `gengetopt`
  is a real Debian/Ubuntu package (`apt-get install gengetopt`) -
  asserted from general knowledge, not verified against an actual apt
  index in this environment (no `apt` available on this macOS dev
  machine) - Phase B should confirm it live in CI rather than trust
  this claim blindly.
- **macOS** (`macos-latest`): no dependencies installed by GitHub by
  default. This project's current default build already finds Apple's
  system libxml2 (`find_package(LibXml2)` resolves to the SDK's
  `.tbd` stub, confirmed locally) and Accelerate for BLAS/LAPACK, with
  zero extra flags - simpler than the deleted workflow's
  Homebrew-OpenBLAS approach. Whether dynamically linking against the
  OS's own libxml2 dylib is portable enough (every "reasonably modern"
  Mac ships it - libxml2 has been part of the OS since early Mac OS X)
  or whether static linking is actually required is an open question
  below, not assumed either way.
- **Windows** (`windows-latest`): the real gap, per the diagnosis
  above. No system libxml2, no system BLAS/LAPACK, no system
  gengetopt. Needs `vcpkg` (or an equivalent) for all three - the
  deleted workflow already had the right *shape* of idea for BLAS
  (`vcpkg` + `CMAKE_TOOLCHAIN_FILE`) but never extended it to
  libxml2/gengetopt, which is the most likely reason it stalled.

## Recommended approach

`vcpkg`, used only where actually needed (Windows, for everything
system-level that isn't otherwise available), with each platform
otherwise preferring its own native/already-proven-working
dependencies - `apt` on Linux (matching the existing job exactly),
Accelerate on macOS (matching the existing local/default build) -
rather than forcing one uniform toolchain (e.g. OpenBLAS everywhere)
across all three OSes, which is what most likely made the deleted
attempt harder than it needed to be.

## Open questions - SETTLED 2026-08-07

- **gengetopt**: resolved, moot (see below, unchanged from before).
- **libxml2 dependency, resolved a different way than expected**: a
  side investigation prompted by "what are the prospects this builds
  right away on someone else's computer?" (`project_xml_optional_
  phylip_default`) found `WITH_LIBXML=ON` used to unconditionally
  compile the XML sources even when `find_package(LibXml2)` found
  nothing - now fixed to auto-detect and gracefully build without XML
  support (a configure-time warning, not a build failure) when
  libxml2 isn't present. Combined with `fastprot`/`fastdist -O` and
  `fnj -I` now defaulting to `phylip` (no external dependency at all)
  instead of `xml`, **libxml2 is no longer a real blocker on any
  platform, including Windows** - a Windows release build can simply
  not provision libxml2 via vcpkg at all and still produce a fully
  functional binary (just without `-O`/`-I xml`), rather than needing
  vcpkg to supply it. This directly simplifies Phase D below.
- **Linking model: dynamic.** Each platform uses its natural/available
  libraries - `apt` on Linux (already proven), Accelerate + optional
  system libxml2 on macOS (already proven), `vcpkg`-provided DLLs
  bundled alongside the `.exe` on Windows for BLAS/LAPACK only (see
  above - not libxml2). No full static linking; `STATIC`'s old intent
  is not being revived.
- **OS version floor**: `*-latest` runners for the push/PR
  build-and-test matrix (Linux/macOS/Windows) - matches the existing
  Linux job, moves forward automatically, and is the right call for a
  job whose purpose is catching breakage against *current*
  environments. The tag-triggered **release** job is the one exception:
  Linux release binaries build on **`ubuntu-22.04`, not
  `ubuntu-latest`** - glibc symbol versioning is forward-only (a binary
  built against a newer glibc won't run on an end-user machine with an
  older one), so the release build should target an older, still-
  supported base rather than whatever's newest on the runner. macOS and
  Windows release builds stay on `*-latest` - macOS's version floor is
  handled via an explicit deployment-target flag instead (see Phase C),
  and Windows doesn't have glibc's forward-only problem in the same way.
- **Trigger and versioning, settled**: two separate jobs, not one.
  Build-and-test matrix (Linux/macOS/Windows) runs on every push/PR,
  same as today's Linux-only job, catching cross-platform breakage
  continuously. A separate release/packaging job runs only on `v*` tags
  - built, not just re-triggered from the push/PR job's artifacts, so
  the Linux leg of *that* job can use `ubuntu-22.04` while the push/PR
  matrix stays on `ubuntu-latest`.
- **Version source, settled (folded in after "is the version noted
  somewhere?", then refined after "is there a way of getting that
  tag-dependent or vice versa?")**: found three *inconsistent* version
  numbers today - root `CMakeLists.txt`'s `PACKAGE_VERSION` (1.0.3,
  CPack-only), `fastprot`/`fastdist`/`fnj`'s independently-hardcoded
  `--version` strings (1.0.10 each, already known/flagged as mismatched
  with `PACKAGE_VERSION` during `gengetopt_migration_plan.md` and
  deliberately left alone then as out of that migration's scope), and
  `fastprot_mpi`'s own hardcoded string (1.0.0, a third value). No
  GitHub release has ever actually been cut (`CHANGELOG.md`'s top
  section is still `## Unreleased`), so nothing downstream depends on
  any of these today - safe to unify now rather than freeze in place.

  **Decision - tag drives version, not the other way round.** The
  first draft ("assert the pushed tag matches a checked-in
  `PACKAGE_VERSION`, fail if not") still left a manual sync step
  someone could get wrong. Better: the release job reads the pushed
  tag (`github.ref_name`, e.g. `v1.2.0`, strip the leading `v`) and
  passes it straight to CMake as `-DPACKAGE_VERSION=1.2.0`, overriding
  whatever's checked in - nothing to keep in sync, no mismatch
  possible because there's nothing left to compare against. For
  everyday (non-release) builds with no tag to derive from - push/PR
  CI, a dev's local `cmake && make` - fall back to `git describe --tags
  --always --dirty` when a `.git` directory is present (informative
  dev-build strings like `1.0.3-42-gabc1234`), falling back further to
  the checked-in static `PACKAGE_VERSION` only when there's no `.git`
  at all (e.g. a `make package_source` tarball, extracted and built
  later with no git history present). The four apps' hardcoded
  `*_VERSION` constants get replaced with the resolved value flowed
  through `config.h` (the same mechanism already carrying
  `WITH_LIBXML` from CMake into the apps' `main.cpp` files - no new
  plumbing pattern needed) - `--version` output is always correct for
  however the binary was actually built, release or not. Mechanical
  fix, not a design question - implemented as part of Phase E (where
  "what version does this release get" naturally belongs), not
  deferred further, but no code changes in Phase A itself.

**New finding while settling these, worth flagging for Phase C**:
GitHub's `macos-latest` runner has been Apple Silicon (arm64) since
the `macos-14` image became its default - meaning the macOS CI job
will need `simde` (this project's non-x86 SSE2-via-NEON shim,
`DNA_b128/sse2_wrapper.h`), unlike what "no additional dependency
installation needed" below originally assumed. `brew install simde`
is one extra line, not a design question, but worth calling out before
Phase C so it isn't discovered as a surprise CI failure instead.

## Phases

### Phase A - DONE (settle the open questions above)

No code changes. See "Open questions - SETTLED" above for the full
record: dynamic linking, `*-latest` for testing / `ubuntu-22.04` for
the Linux release leg, two separate workflows, tag-drives-version.

### Phase B - moot

Was going to add `apt-get install gengetopt` to remove the FTP
self-build dependency from every CI run - no longer needed,
`gengetopt_migration_plan.md` removed gengetopt from the build
entirely, so this CI job had nothing left to fix here.

### Phase C - DONE (add macOS to the matrix)

`macos-latest` added to `build-and-test.yml`'s matrix. Only new
dependency needed was `simde` (`brew install`) - Accelerate and system
libxml2 already sufficient, matching the dynamic-linking decision.
Verified green immediately (this session's own dev machine is Apple
Silicon, so this path had effectively already been exercised all
along). Also switched the build step to `cmake --build build
--parallel` (portable) instead of `-j"$(nproc)"` (Linux/GNU-only,
would have failed on macOS and Windows both).

### Phase D - DONE (add Windows to the matrix)

The real work, and it was real: 12 CI iterations to reach fully green,
with no local Windows machine to test any of it against beforehand.
`vcpkg` in manifest mode (`vcpkg.json` - `lapack-reference` only, not
libxml2, per Phase A) via a relative `CMAKE_TOOLCHAIN_FILE` path
(`${{ github.workspace }}` is not reliably available inside
`strategy.matrix` - found the hard way on attempt 1), Ninja instead of
the default Visual Studio generator (single-config output path,
matching Linux/macOS), `ilammy/msvc-dev-cmd` for the compiler
environment, and vcpkg binary caching via the GitHub Actions cache
(`x-gha`) so `lapack-reference`'s ~15-minute cold Fortran build only
happens once, not on every run.

Six distinct, real source-level MSVC/Windows portability bugs found
and fixed along the way - none of them CI plumbing, all genuine
first-time-ever-Windows-compiled defects:

1. `BitVector.cpp`: `unsigned long &` bound to a `vector<size_t>`
   element - silently identical types on Linux/macOS (LP64), genuinely
   different types on Windows (LLP64). Fixed with `auto &`.
2. `ambiguity_nucleotide.hpp`/`DNA_b128_String.cpp`: 22 occurrences of
   C99 compound-literal syntax (`(Type){...}`), never valid standard
   C++, silently accepted by GCC/Clang as an extension. Fixed by
   dropping the invalid parens (`Type{...}`, standard brace-init).
   Also removed `sse2_wrapper.h`'s dead `getticks()` (GCC inline
   assembly, zero callers, no MSVC equivalent).
3. `NeighborJoining.hpp`: 10 C99 Variable-Length-Array declarations
   sized by a runtime value, across 6 function templates - never valid
   standard C++ either, same GCC/Clang-extension story. Fixed with
   `std::vector<T>` - the more idiomatic modern C++ answer anyway, not
   just a portability patch.
4. `ProtSeqCompare.cpp`: `__builtin_popcount` (GCC/Clang builtin, no
   MSVC equivalent) in genuinely hot-loop SIMD comparison code - fixed
   with a small `popcount16()` wrapper (`__popcnt` from `<intrin.h>` on
   real MSVC, the GCC builtin everywhere else including clang-cl).
5. `rand_r()` (POSIX-only) in `ProtSeqCompare_test.cpp` and
   `bench_primitives.cpp` - switched to `std::mt19937`, standard C++11
   and more reproducible than `rand_r()` ever was.
6. `isatty()`/`STDIN_FILENO` (POSIX-only, `<unistd.h>` does not exist
   on MSVC) in `fnj`/`fastprot`/`fastprot_mpi`'s "no input data" stdin
   check - fixed with a small `stdinIsATerminal()` helper using MSVC's
   `_isatty()`/`_fileno()` from `<io.h>`.

Plus real build/link/runtime-behavior gaps, not just compile errors:

7. `fnj/DataInputStream.hpp` unconditionally included `<libxml/tree.h>`
   with zero actual use of any libxml type - a real bug already, just
   never observed before because libxml2 was always present everywhere
   this project had ever been built. Deleted (not guarded) since it
   was dead weight regardless of platform.
8. `target_link_libraries(... m ...)`: Windows has no separate math
   library at all (`m.lib` does not exist - every math function is
   already part of the CRT). Introduced `FASTPHYLO_MATH_LIB`, `m` only
   under `if(UNIX)`.
9. Windows' CRT defaults stdout (and files opened via `fopen(..., "w")`/
   `ofstream` in text mode) to text mode, silently translating every
   `\n` this project's code writes into `\r\n` - invisible in the
   compile/link/unit-test steps, but broke every single one of
   `RunExamples.sh`'s byte-exact comparisons at once. Fixed at both
   ends: `file_utils.cpp`'s four write-opening helpers now open in
   binary mode, and a `setStdoutBinaryMode()` helper (`_setmode()`)
   runs first thing in each app's `main()`.
10. Two examples (`ex3`/`ex9`) need real XML support, which this
    Windows build genuinely does not have (Phase A's own deliberate
    decision); one (`ex16`, a raw binary float dump) is not required by
    IEEE 754 to be bit-identical across math library implementations
    for transcendental functions like `log`/`exp`. Both are real,
    inherent, permanent platform differences, not bugs - excluded from
    the byte-exact check specifically on the Windows leg via
    `RunExamples.sh`'s new `SKIP_EXAMPLES` env var (a no-op everywhere
    else), verification staying full-strength (all 20 examples) on
    Linux/macOS.

Also caught and fixed a bug in the fix for #10 itself: the first
attempt added a second, Windows-only diff loop *after* `bash
RunExamples.sh` in the workflow step - but `RunExamples.sh` already had
its own internal, unconditional diff-and-report loop at the very end
(read early in this engagement, not accounted for when writing that
fix), which ran first and failed the whole step before the new loop
was ever reached. The real CI run's output (unchanged from before the
supposed fix) is what caught it. Fixed properly by parameterizing
`RunExamples.sh`'s own existing loop instead of duplicating it.

Confirmed green with `continue-on-error` removed (it was only there
while this phase was being stabilized).

### Phase E - Packaging and release automation

Only once B-D all build and pass `RunExamples.sh` reliably on every
push/PR, not just at tag time: per-OS artifact packaging (zip/tar.gz),
`actions/upload-artifact` for regular builds, `softprops/action-gh-release`
(or equivalent) on tagged releases - matching the deleted workflow's
shape, now with working dependency provisioning underneath it.

### Phase F - Final verification

A real tagged test release, confirming the produced binaries actually
run standalone per OS - particularly checking that nothing
dynamically-linked on the Windows build silently depends on the build
runner's own `vcpkg`/toolchain paths being present at runtime.

## Decisions

- Builds on the existing, working `build-and-test.yml` rather than
  reviving the deleted `release.yml` wholesale - same lesson as the
  `STATIC` removal (`lint_plan.md`'s Phase 6): fix the actual gap
  (dependency provisioning, most acutely on Windows) rather than debug
  an abandoned attempt's accumulated complexity from where it left off.
- `STATIC` (the old self-built libxml2/zlib/gengetopt toolchain,
  removed in `lint_plan.md`'s Phase 6) is *not* the mechanism this plan
  revives - `vcpkg` (or platform-native package managers) replaces its
  role, with current dependency versions and no FTP reliance.
- Branch/timing: not blocking `lint_plan.md`'s remaining Phase 7, and
  not started - this doc exists to capture the investigation and
  decision to defer, per direct discussion.
