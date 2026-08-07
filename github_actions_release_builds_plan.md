# Plan: GitHub Actions cross-platform release builds

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

### Phase A - Settle the open questions above

No code changes. Answers here shape every later phase - in particular
whether Windows/macOS need full static linking (closer to reviving
`STATIC`'s original *intent*, properly this time) or whether dynamic
linking against each OS's own libraries is good enough.

### Phase B - Fix `build-and-test.yml`'s Linux job in place (moot)

Was going to add `apt-get install gengetopt` to remove the FTP
self-build dependency from every CI run - no longer needed,
`gengetopt_migration_plan.md` removed gengetopt from the build
entirely, so this CI job has nothing left to fix here.

### Phase C - Add macOS to the matrix

`macos-latest`, no additional dependency installation needed if Phase
A settles on dynamic linking (Accelerate + system libxml2 already
sufficient, per the current default build) - same
configure/build/ctest/RunExamples.sh shape as the Linux job.

### Phase D - Add Windows to the matrix

The real work: `vcpkg` in manifest mode (a `vcpkg.json` listing
libxml2 and a BLAS/LAPACK provider - gengetopt no longer needed here,
removed entirely by `gengetopt_migration_plan.md`), `CMAKE_TOOLCHAIN_FILE`
wiring, and verifying `RunExamples.sh` actually runs correctly under
Windows path/line-ending conventions - this project's file I/O has had
real CRLF-handling surprises before (see project memory from the lint
pass), worth specifically checking here rather than assuming it's fine.

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
