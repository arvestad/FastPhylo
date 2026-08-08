# Plan: catalogue of `core/` modernization debt

## Goal

Direct instruction: this needs to end with genuinely *modern* C++, not
just C++ that compiles cleanly - and a catalogue first, decisions and
execution second, rather than chasing each finding ad hoc as it turns
up. Grew out of `legacy_error_handling_plan.md` (all 5 phases done)
once it became clear `log_utils.hpp`, `InitAndPrintOn_utils.*`,
`Object.*`, and `stl_utils.*` all have the same character: real, live
code sitting next to a large amount of dead code, written in a
pre-C++11, sometimes Java-flavored style throughout.

Everything below was verified directly (grep for real call sites,
cross-checked against what's actually built) on 2026-08-08/09, not
assumed from reading alone.

## "Is `FloatDistanceMatrix.cpp` needed?" - answered, and acted on (2026-08-08/09)

**No, and neither was anything else in its parallel "float precision"
NJ/tree pipeline - all deleted, `fa9aebb`/`b36b725`.** Confirmed by
tracing the actual call graph, not just grepping the type name:
`fnj/main.cpp` only ever instantiates the `buildTrees<T>` template
with `T=StrDblMatrix`, so `computeNJTree(StrFloMatrix&, ...)` was
never reached; `DataInputStream`'s abstract interface no longer even
declares a `StrFloMatrix` `readDM` overload (commented out, with a
2016-06-14 note from Lasse explaining why); `BinaryInputStream.cpp`
already had its sibling `StrFloMatrix` overload deleted in this
session's earlier `fnj_binary_input_gap` work for the identical
reason ("had zero callers anywhere in the codebase"); every
`SequenceTree::*Float*` method and every `computeFloat{NeighborJoining,
BioNJ,FNJ}Tree` in `NeighborJoining.hpp` had zero external callers.
Deleted ~1200 lines total: `FloatDistanceMatrix.hpp`/`.cpp`/`_impl.hpp`
outright, the `computeFloat*Tree` functions, all `SequenceTree::*Float*`
methods, and the one remaining unreachable stub
(`XmlInputStream::readDM(StrFloMatrix&, ...)`, which just printed "Not
implemented!" and called `exit(-1)` - not part of the virtual
interface, so nothing could ever call it). `StrFloRow`
(`DistanceRow<...,float,...>`, a completely separate class hierarchy
used by `fastdist`'s real, live memory-efficient row-streaming output)
is unrelated and was untouched.

## `Extrainfos.hpp` - not dead, but a real interface-design smell

`using Extrainfos = std::vector<std::string>;` - a bare, semantically
opaque type alias threaded through the `readDM`/`printStartRun`
signatures of **every** `DataInputStream`/`DataOutputStream` subclass
across all three formats (Phylip, binary, XML) - roughly 60 files
reference it. But only **XML** ever populates or reads it: it's
per-sequence-name-indexed freeform text serialized into `<identity>`
tags (`XmlOutputStream.cpp`) and parsed back (`XmlInputStream.cpp`).
Phylip and binary formats carry an always-empty `Extrainfos&` parameter
through their entire call chain without ever touching it.

**Verdict: needs-design-decision, not mechanical.** The problem isn't
size (the feature itself is small) - it's that a single XML-specific
concern was baked into the generic cross-format interface every
`DataInputStream`/`DataOutputStream` subclass has to implement, plus a
name that gives no hint what it actually holds. A real fix means
either renaming it to something that describes what it is (e.g.
`SequenceIdentityAnnotations`) and/or removing it from the generic
interface entirely (e.g. an XML-only side channel rather than a
parameter every format's `readDM`/`printStartRun` must carry). Not
attempted here - flagged for the same future design-decision pass as
`Object.hpp`, given both require touching a wide interface rather than
a local, mechanical change.

## "What is the purpose of `SequenceTree.cpp`?" - answered first

**It's genuine, central, actively-used code - not cruft.** It's the
tree data structure `fnj` (Fast Neighbor Joining) builds and returns:
a specialization of the generic `Tree<Data, DataInit, DataPrintOn>`
template (`Tree.hpp`) with `Sequence_double` (a sequence + an edge
length) as the per-node payload. Confirmed live via direct callers:
`fnj/main.cpp`, `fnj/TreeTextOutputStream.hpp`, `fnj/DataOutputStream.hpp`,
`NeighborJoining.cpp`/`.hpp` (the actual NJ algorithm), `Sequences2DistanceMatrix.hpp`,
`phylip_interleaved_reader.hpp`, and directly exercised by
`tests/SequenceTree_PhylipReader_test.cpp` (one of the 4 real `ctest`
suites). It implements real phylogenetics operations: Robinson-Foulds
distance between two trees, tree-to-distance-matrix conversion,
log-likelihood computation, canonicalization for tree comparison,
degree-2 node contraction, sequence-to-node mapping. None of that is
in question here - the algorithms stay untouched.

It does have its own, narrower modernization debt (catalogued below
with everything else): C-style field-access macros, and one dead
overload.

## Catalogue

Each item: **verdict** - dead (delete outright), keep-but-modernize
(real code, wrong-era style), or needs-design-decision (too large/
cross-cutting for a mechanical pass).

### `log_utils.hpp` - a second, uncoordinated error/logging system

**Done, 2026-08-09.** Full file already read together earlier. 11
macros total.

- **Dead (delete outright)**: `PRINT_V`, `PRINT_TIME`, `PRINT_EXP`,
  `LINE`, `SEPARATOR` - confirmed 0 real call sites, verified
  excluding comments. **Correction to this catalogue's earlier
  claim**: `PRINT` was *not* actually dead - one real, non-commented
  call site survived in `Tree_impl.hpp:715`
  (`PRINT(n->parent == NULL)`, inside an `if (n->parent == NULL)`
  block - printing a condition already established by the enclosing
  `if`, zero debugging value). The original "0 real call sites" claim
  was a genuine miss in the audit, only caught while executing this
  phase (not before deleting anything - the macro was still fully
  defined at that point, so nothing broke). Deleted that one call site
  outright rather than keep the whole macro family alive for it.
  (`ASSERT`/`ASSERT_OP`/`DEBUG` were already deleted in
  `legacy_error_handling_plan.md`'s Phase 5 - one of them, `DEBUG`,
  had a latent syntax bug that only compiled because it was never
  invoked; same pattern here.)
- **`MEM_CHECK` (2 call sites)** - inlined a direct `if (ptr ==
  nullptr) { std::cerr << "out of memory" << std::endl; std::abort();
  }` at both `DNA_b128_String.cpp` call sites, matching the reasoning
  above (minimize further allocation while already handling OOM).
  Macro deleted.
- **`USER_ERROR`, `PROG_ERROR` (confirmed via full-tree re-grep before
  touching anything: 13 + 12 = 25 real call sites across 12 files)** -
  every call site renamed to `THROW_EXCEPTION` (syntax matched
  exactly, purely mechanical). Both macros deleted.
- **`USER_WARNING` (14 real call sites across 6 files)** - every call
  site replaced with a direct `std::cerr << "warning: " << ... <<
  std::endl;` - no function/macro wrapper at all, since (unlike
  `THROW_EXCEPTION`) it needs no control-flow beyond an ordinary
  statement. Dropped the old boxed format's file/line/function context
  deliberately - these are user-facing messages about *data*
  (non-finite floats, mismatched tree leaf sets), not internal
  diagnostics, so that context wasn't earning its keep for the actual
  audience. Macro deleted.
- **Scope**: 15 files use `USER_ERROR`/`PROG_ERROR`/`USER_WARNING`/
  `MEM_CHECK` combined - `core/` (`DistanceMatrix.cpp`, `Sequence.cpp`,
  `SequenceTree.cpp`, `InitAndPrintOn_utils.hpp`, `stl_utils.hpp`,
  `Tree.hpp`), `dna/` (`DNA_b128_String.cpp`, `NeighborJoining.cpp`,
  `Sequences2DistanceMatrix.cpp`, `string_compare.cpp`,
  `ambiguity_nucleotide.hpp`, `NeighborJoining.hpp`), `io/`
  (`PhylipDmOutputStream.cpp`, `XmlOutputStream.cpp`), and one
  `fastprot_mpi` file (`PhylipDmOutputStream.cpp` - unverified, no MPI
  available, same standing caveat as everywhere else in this
  engagement).

### `InitAndPrintOn_utils.cpp`/`.hpp` - "shared utils" that mostly isn't shared

- **Dead (delete outright)**: `string_int`, `string_double`,
  `int_double` structs and their `operator>>`/`operator<<` pairs - 0
  real callers anywhere, roughly 90 of the `.cpp`'s 137 lines.
- **Keep, but relocate**: `Sequence_double` + its `operator>>`/
  `operator<<` - live, but its only consumers are `NeighborJoining.cpp`/
  `.hpp` and `SequenceTree.hpp`. Not a shared utility in any real
  sense - one consumer cluster wearing a generic-sounding filename.
  **Verdict**: move it to live next to `SequenceTree`/`NeighborJoining`
  rather than a file named for a category that's 3/4 fictional.
- **Keep as-is, genuinely shared**: `Data_init<T>`/`empty_Data_init<T>`/
  `Data_printOn<T>`/`empty_Data_printOn<T>` template functors - used
  across `Tree.hpp`, `DistanceMatrix.hpp`, `DistanceRow.hpp`,
  `FloatDistanceMatrix.hpp`, `SequenceTree.hpp`, `NeighborJoining.hpp`.
  This part earns the "utils" framing; unrelated to the dead parsing
  structs above.

### `Object.hpp`/`.cpp` - not modern, but not dead - needs a design decision

A hand-rolled, Java-style common base (`virtual std::string toString()`,
`virtual size_t hashCode()`, `virtual bool equals(const Object*)`,
`virtual operator==`), all dispatched through a vtable. Modern C++
doesn't need a common base class for any of this: `operator<<` is
found via ADL per-type with no base class involved, equality is
`operator==` defined per-type, hash-map keys use `std::hash<T>`
specializations rather than a virtual `hashCode()`. Structurally, yes -
this predates (and costs more, via virtual dispatch, than) how it'd be
written today.

**But it's not dead code** - confirmed live: `equals()`/`hashCode()`
via the `objeq`/`objhash` functor structs (`SequenceTree`'s hash
containers), `objInitFromStream()` overridden and called across 6
classes (`Sequence`, `DistanceRow`, `DistanceMatrix`,
`FloatDistanceMatrix`, `Tree`, `BitVector`), and `printOn()`/the free
`operator<<(ostream&, const Object&)` used throughout for diagnostic
output (including, until recently, `Exception` - see
`legacy_error_handling_plan.md`).

**Verdict: needs-design-decision, not a mechanical pass.** A real
modernization means redesigning how 7 classes handle printing,
equality, hashing, and stream-init simultaneously - touching every
piece of code that builds a hashmap/hashset keyed by one of them. Too
large and cross-cutting to bundle into this catalogue's mechanical
items; flagged for a dedicated future decision on its own, once the
smaller items below are done and the remaining shape of the problem is
clearer.

### `stl_utils.hpp`/`.cpp` - several more dead operators, one redundant-but-live one

- **Likely dead, needs a build-time check to confirm** (operator-
  overload call sites don't grep reliably by text search alone - no
  evidence found via search, but that's not proof of absence the way
  it is for named functions): 4 `operator+(string, int/float)`
  overloads, built via fixed-size `snprintf` buffers with leftover
  `#ifdef WIN32`/`_snprintf` handling. Even if a live call site turns
  up: this is exactly the kind of implicit, surprising operator
  overloading modern style avoids (`"x" + 5` silently means "stringify
  and concatenate") - `std::to_string` is the direct modern
  replacement regardless.
- **Likely dead, same caveat**: `operator>>`/`operator<<(std::string*)`
  - both just delegate to the standard reference-based operators; no
  evidence of any caller streaming through a raw `string*` rather than
  a reference.
- **Likely dead, same caveat, plus real bugs even if kept**:
  `operator>>`/`operator<<` templates for `std::vector<T>`. No
  evidence of any caller. Independent of dead-or-not: declared in the
  *global* namespace for a *standard-library* template - risks silent
  conflicts with any other header doing the same; the `<<` overload
  takes `vector<T>&` where it should be `const&`; `vec[vec.size()-1]`
  underflows unguarded if `vec` is empty.
- **Dead together**: `ltstr` + `str2int_map` + `print_map` - `ltstr`
  has 0 callers, and the other two only exist to use it.
- **Live but redundant**: `hashstr` (used by `Sequence.cpp`/
  `SequenceTree.cpp`) duplicates what the standard library already
  gives `std::unordered_map<std::string, V>` for free - no custom
  hash functor is needed for plain `std::string` keys at all.
  `eqstr` similarly duplicates `std::equal_to<std::string>`
  (default). **Verdict**: once confirmed no longer needed by any
  hashmap declaration, delete; `str2int_hashmap`/`str2str_hashmap`
  (only the former is live) can drop to plain
  `std::unordered_map<std::string, int>` with no custom functors.
- **Dead**: `obj_ptr2obj_ptr_hashmap`/`obj2obj_hashmap` - their only
  reference anywhere was `SequenceBasedNJ.cpp`, itself confirmed
  unbuilt dead weight (not registered in any `CMakeLists.txt`).
- **Live, genuinely fine**: `appendUntil` (the 4-arg raw-buffer
  version - not to be confused with the dead 3-arg `istream`-based
  `appendUntil` already deleted from `file_utils.hpp` in
  `legacy_error_handling_plan.md`'s Phase 4) and `appendAllNonChars` -
  both used by `Sequence.cpp`'s real FASTA-parsing code. Not touched
  here; C-style buffer parsing is defensible in a hot parsing loop,
  and neither is currently broken.

### `SequenceTree.cpp`/`.hpp` - core algorithms stay, two narrow items

- **Dead**: `mapSequencesOntoTree(char **nameseqPairs, int numPairs)` -
  0 callers anywhere; only the `istream&` overload is used. Delete the
  raw-C-array overload.
- **Dead**: the `SEQUENCE(node)` access macro - 0 occurrences outside
  its own `#define`.
- **Keep, but modernize**: `NAME(node)`/`SEQ(node)`/`EDGE(node)` -
  97 combined call sites (31 + 23 + 43) of C-preprocessor field-access
  macros standing in for what should be member accessor methods
  (`node->name()`/`node->seq()`/`node->edgeLength()` or similar). Real
  code, real usage, wrong era - a big, mostly-mechanical rename/
  refactor given the call-site count, not a redesign.

## Suggested execution order

1. **Done.** Dead-code deletions across every file above: the
   `log_utils.hpp` macro family, `string_int`/`string_double`/
   `int_double`, `ltstr`/`str2int_map`/`print_map`,
   `obj_ptr2obj_ptr_hashmap`/`obj2obj_hashmap`, `str2str_hashmap`,
   `mapSequencesOntoTree(char**, int)`, `SEQUENCE(node)`. The 4
   `operator+` overloads and the `string*`/`vector<T>` stream
   operators were confirmed dead (all but one) via an actual
   build-time removal check (`#if 0` + full rebuild), not just grep -
   this surfaced the one genuinely live `operator+(string, int)`
   overload (2 call sites), which was replaced with `std::to_string`
   instead of kept.
2. **Done.** `USER_ERROR`/`PROG_ERROR` -> `THROW_EXCEPTION`,
   `USER_WARNING` -> direct `std::cerr`, `MEM_CHECK` inlined -
   mechanical, touched 16 files. Verified with a full rebuild +
   `ctest` + `RunExamples.sh` + `RunCliChecks.sh` on both Clang and a
   real local GCC 14 build (per the `<iomanip>` lesson - these are all
   widely-included headers), plus `-DWITH_LIBXML=OFF`.
3. **Done.** `Sequence_double` stays in `InitAndPrintOn_utils.*` for
   now (not relocated) - `string_int`/`string_double`/`int_double`
   were deleted as part of step 1, which was the main win; relocating
   the one remaining live struct is cosmetic and lower priority than
   it looked before the dead 90% was cut.
4. **`NAME`/`SEQ`/`EDGE` macros -> accessor methods** - not started.
   Bigger (94 combined call sites now that `SEQUENCE`'s dead uses are
   gone) but still mechanical; the natural next step.
5. **`Object` redesign** - explicitly out of this catalogue's
   execution; needs its own dedicated discussion given its size and
   how many classes depend on it (`Sequence`, `DistanceMatrix`,
   `BitVector`, `DistanceRow`, `FloatDistanceMatrix`, `Tree`,
   `TreeNode`, plus every hashmap keyed by one of them).

## Explicitly out of scope

- `Object.hpp`/`.cpp`'s full redesign (see above) - catalogued, not
  scheduled.
- `fastprot_mpi` - unverified throughout (no MPI available); changes
  there follow by inspection/analogy only, same standing rule as the
  rest of this engagement.
- Any change to `Tree.hpp`/`Tree_impl.hpp` beyond the `NAME`/`SEQ`/
  `EDGE` macro accessors - those files are large (464 + 1627 lines)
  and haven't been audited line-by-line the way `core/`'s smaller
  files have; a full pass there would be its own catalogue.
