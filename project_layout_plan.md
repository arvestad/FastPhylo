# Plan: Restructure FastPhylo to a library + apps layout

## Goal

Give the project a real library boundary: a single, discoverable public
`include/` directory, implementation that mirrors it, and per-program
(`apps/`) directories that are thin consumers instead of near-independent
copies of shared infrastructure. Same discipline as `plan.md`
(speed2026a) and `modernization_plan.md`: verify byte-identical output
(`RunExamples.sh`, `ctest`) after every step, never trade correctness
for structure.

Branch: TBD (off `modernize-cpp17`'s tip, so it carries forward all
prior work rather than starting from `master`).

## Why

Lasse's own words: "I have had a hard time navigating the project."
Two concrete, evidence-backed problems, not just aesthetics:

1. **No discoverable public API.** Headers and implementation are flat
   and interleaved across `src/c++/`, and every program's include path
   is opened wide to all of it (`target_include_directories(fastphylo
   PUBLIC ...)` from this plan's own Phase 5 has to list five separate
   directories). There's no single place to look and see "this is what
   the library exposes."
2. **The same infrastructure is copy-pasted 3-4 times, and it rots
   independently.** `programs/fastdist/`, `programs/fnj/`,
   `programs/fastprot/`, and `programs/fastprot_mpi/` each carry their
   own near-identical copies of `DataInputStream`/`DataOutputStream`,
   `Extrainfos.hpp`, and the Xml/Phylip/Binary/Fasta stream readers and
   writers. This is not hypothetical risk - it's the direct cause of a
   bug pattern already found and fixed **three separate times** this
   engagement: `fastprot`'s `BinaryDmOutputStream` had a heap-corruption
   bug (`delete` on a non-`new`'d pointer), and `fastdist`'s
   independently-written copy of the *same class* had a *different* bug
   (a file-descriptor leak) from the *same* underlying confusion about
   `fp`'s ownership. Two people (or two moments) made two different
   mistakes implementing the same idea twice, because it was implemented
   twice. `fastprot_mpi` additionally appears to have fallen behind
   `fastprot` itself - it has no copy of `ProtSeqCode.cpp`/
   `ProtSeqCompare.cpp` at all, meaning it never picked up speed2026a's
   SIMD work and is running the old scalar path silently.

## Proposed layout

```
FastPhylo/
├── CMakeLists.txt
├── include/
│   └── fastphylo/              # the ENTIRE public API - nothing else
│       ├── core/
│       ├── dna/
│       ├── io/
│       └── protein/
├── src/
│   └── fastphylo/               # implementation, mirrors include/fastphylo/ 1:1
│       ├── core/  dna/  io/  protein/
├── apps/                         # renamed from programs/; thin consumers only
│   ├── fastdist/    main.cpp + gengetopt/ + a short CMakeLists.txt
│   ├── fnj/
│   ├── fastprot/
│   └── fastprot_mpi/
├── tests/                        # pulled out of src/c++/code_tests/
├── benchmarks/
├── examples/
└── docs/                         # renamed from src/docbook/
```

## Concrete file mapping (audited against the current tree, not guessed)

**`core/`** (used by both the DNA and protein sides; no domain-specific
dependencies) - `BitVector`, `Exception`, `InitAndPrintOn_utils`,
`Object`, `Sequence`, `SequenceTree`(`_MostParsimonious`), `Tree`
(`_impl`), `Simulator`, `arg_utils_ext`, `arg_utils` (C), `file_utils`,
`stl_utils`, `std_c_utils` (C), `xml_output_global`, `DistanceMatrix`
(`_impl`), `FloatDistanceMatrix` (`_impl`), `DistanceRow` (`_impl`),
`log_utils.hpp`, `nucleotide.hpp`. (`nucleotide.hpp` looks DNA-specific
by name, but it's actually included by `Sequence.cpp`/`SequenceTree.cpp`
/`SequenceTree_MostParsimonious.cpp` - core files, not just the DNA
side - confirmed by grep, so it belongs here, not in `dna/`.)

**`dna/`** - `DNA_b128/` (4 files), `distance_methods/`
(`LeastSquaresFit`, `NeighborJoining`), `sequence_likelihood/`
(`Kimura2parameter`, `TamuraNei`, `ambiguity_nucleotide`,
`dna_pairwise_sequence_likelihood`, `string_compare`).

**`io/`** - `DataInputStream`, `DataOutputStream`, `Extrainfos`,
`XmlInputStream`/`XmlOutputStream`, `PhylipMaInputStream`/
`PhylipDmInputStream`, `PhylipDmOutputStream`, `FastaInputStream`,
`BinaryDmOutputStream`/`BinaryInputStream`, `TreeTextOutputStream`,
`fileFormatSchema.hpp`. **This is not a file move** - it's a
consolidation. Each of these currently exists as 2-4 near-but-not-quite
identical copies across `fastdist`/`fnj`/`fastprot`/`fastprot_mpi`
(different because they read different things: sequences vs. distance
matrices vs. trees). Every copy needs to be diffed against its siblings
before merging, so real divergence (not just formatting drift) survives
as an explicit specialization instead of being silently averaged away.

**`protein/`** - `ProtSeqCode`, `ProtSeqCompare`, `Matrix`,
`MaximumLikelihood`, `ExpectedDistance`, `ProtDistCalc`, `ModelMatrix`,
`ProtSeqUtils`, sourced from `programs/fastprot/` only. `fastprot_mpi`
is excluded from this plan entirely for now (see "Deferred:
`fastprot_mpi`" below) - it keeps using its own existing, untouched
copies of these files until that work happens.

**`apps/`** - `fastdist`, `fnj`, `fastprot` (not `fastprot_mpi` - see
below): each keeps `main.cpp` and its `gengetopt/` directory. Once
`io/` and `protein/` are consolidated into the library, whatever's left
in each app directory should genuinely just be CLI wiring plus true
per-program logic - if a file doesn't shrink to near-nothing after the
library absorbs the shared pieces, that's a signal it wasn't actually a
duplicate and belongs where it is.

## Phases

Staged so no single commit is both large and risky at once - the
mechanical moves (low risk) happen before the consolidation work (real
engineering, needs case-by-case review like Phase 2's RAII work did).

### Phase A - Scaffolding
Create `include/fastphylo/{core,dna,io,protein}/` and
`src/fastphylo/{core,dna,io,protein}/` empty, wire up
`target_include_directories`/`target_sources` for the new paths
alongside the old ones (both work simultaneously), no files moved yet.
Verify the build is unaffected.

### Phase B - Move `core/` and `dna/`
Pure relocation - both are already single-copy. Move files,
`#include "X.hpp"` -> `#include <fastphylo/core/X.hpp>` (or `dna/`)
everywhere they're used, delete the old empty directories, remove the
old include paths from CMake. Highest-confidence phase - proves the
new layout end-to-end before the hard parts.

### Phase C - Consolidate `io/`
The hard one. Diff each duplicated class across `fastdist`/`fnj`/
`fastprot` (not `fastprot_mpi` - deferred, see below), identify genuine
behavioral differences vs. copy-paste drift, design the shared version
(plain shared implementation where the copies are identical; a small
extension point - virtual method, template parameter, or strategy
object - where they're not). Move the result into `io/`, update the
three apps to use it, delete their duplicates. This is where the
fastprot/fastdist `BinaryDmOutputStream` bug pattern gets structurally
prevented, not just individually fixed again.

### Phase D - `programs/` → `apps/`
Rename, confirm each app directory is now genuinely thin, short
per-app `CMakeLists.txt` files if that reads better than the current
one-big-file approach (open question, not decided yet).

### Phase E - `tests/` and `docs/`
Move `src/c++/code_tests/` to top-level `tests/` (sibling to
`examples/`/`benchmarks/`, where you'd actually look for it); rename
`src/docbook/` to `docs/`. Low risk, mechanical.

### Phase F - Verification and sign-off
Full rebuild + `ctest` + `RunExamples.sh` (paths inside it reference
`FASTPHYLOPATH=.../build/src/c++` - needs updating for the new binary
location) + a check of `.github/workflows/build-and-test.yml` and
`benchmarks/run_benchmarks.sh` for any now-stale paths. Update
`README`/memory notes to match the new layout.

## Decisions (settled 2026-08-01)

- **Bucket names** (`core`/`dna`/`io`/`protein`): confirmed as-is.
- **Namespace**: not doing it now. Main benefit is collision safety if
  the library gets linked into another project (`Object`/`Exception`/
  `Tree`/`Sequence` are generic, collision-prone names) - genuinely
  relevant given the "provided to whatever other project might be
  added" goal, but speculative until there's an actual second consumer,
  and the project is expected to stay well-defined rather than grow
  into something large. `include/fastphylo/` already gives "is this
  ours" clarity in headers without touching ~30K lines of
  implementation. Staying cheap to add later - the directory convention
  is already namespace-shaped.
- **Sequencing vs. linting**: restructuring and linting are orthogonal;
  linting waits until this is done, so lint fixes land in final file
  locations rather than files about to move.

## Deferred: `fastprot_mpi`

Explicitly waiting - "can wait a bit," to be picked up on request
separately from this plan. When it happens, the goal is **folding
`fastprot_mpi` into `fastprot`** (one target, MPI-enabled as a build
option - in the spirit of the existing `WITH_LIBXML`/`BUILD_WITH_MPI`
toggles) rather than giving it its own parallel `protein/`-consuming
app forever. That reframes what was Phase D above: instead of
consolidating shared code between two permanent targets, the end state
has at most one `fastprot` app directory. Also explains the apparent
SIMD staleness noted above (no `ProtSeqCode`/`ProtSeqCompare` copy) -
folding resolves that gap as a side effect rather than needing a
separate reconciliation step.

## What's explicitly out of scope for this plan

- `fastprot_mpi` (see "Deferred" above - separate future request).
- Migrating off `gengetopt` (Lasse asked about this separately -
  CLI11 is the recommendation, but it's an orthogonal CLI-parsing-library
  swap, not a file-layout change; `apps/*/gengetopt/` would naturally
  disappear if that migration happens, but this plan doesn't assume it
  will).
- The linting/code-smell pass - deliberately sequenced after this plan
  (see Decisions above).
- Everything already listed as deferred in `modernization_plan.md`
  (STATIC's FTP-hosted deps, docbook's deprecated CMake API, the
  `Matrix`/Eigen investigation).
