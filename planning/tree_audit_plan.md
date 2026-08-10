# Plan: `Tree.hpp`/`Tree_impl.hpp` audit

## Goal

`Tree<Data,DataInitializer,DataPrintOn>`/`TreeNode<...>` (`Tree.hpp`,
584 lines; `Tree_impl.hpp`, 1621 lines) are the one remaining
`core/` file pair never fully audited this engagement -
`modern_cpp_core_cleanup_plan.md` explicitly carved them out ("large
... and haven't been audited line-by-line the way `core/`'s smaller
files have"), and the `Object.hpp` removal plan only touched the
printing/equality/`objInitFromStream` interface, not the rest of
these files' ~2200 lines.

**Deliberately scoped for after the 2.0beta release, not before.**
Scoping this (2026-08-10) surfaced three narrow, real bugs
(`TREE_ASSERT` running unconditionally in release builds,
`operator=`'s missing self-assignment guard, `makeCanonical`'s
variable-shadowing null-pointer bug) - those were small and contained
enough to fix immediately (commit `f55092f`) rather than wait. The
rest of this plan - a ~14-method dead-code cluster, several structural
cleanups, one error-handling inconsistency - is real but larger and
lower-risk-to-defer; a beta period is exactly the right time to let
more real usage surface anything this survey missed before removing
that much interconnected code at once.

Confirmed via this engagement's own earlier work: `SequenceTree`
(`Tree<Sequence_double, Data_init<Sequence_double>,
Data_printOn<Sequence_double>>`) is the only live instantiation of
`Tree<>` anywhere in the built tree - the only other one,
`b128Tree` in `SequenceBasedNJ.cpp`, has been confirmed unbuilt
(no `CMakeLists.txt` reference) since `legacy_error_handling_plan.md`'s
Phase 0. Every dead-code finding below can be checked against
`SequenceTree`'s actual call graph without worrying about a second,
differently-shaped consumer.

## Finding 1 - a genuinely dead ~14-method cluster, closed on itself

Surveyed (via a dedicated read-through, not just grep) every public
method on both classes for external callers (i.e. callers outside
`Tree.hpp`/`Tree_impl.hpp`/`SequenceTree.hpp`/`SequenceTree.cpp`
themselves). Found a cluster with **zero callers anywhere in the
repo**, not just outside the base-infrastructure files - and the
cluster is closed: every method in it is only ever called by another
member of the same cluster.

- `Tree::isBinary`, `reRootAtHighDegreeNode`, `moveNode`, `splittTree`,
  `joinTreeAtRoot`/`joinTreeAtNode`/`joinTreeAtNodes` (a mutually-
  calling trio with no outside entry point), `insertNodeOnPathToParent`,
  `shortcutNode`, `removeAndDelete` (used only by `shortcutNode`),
  `collapse`/`collapseChildren`, `addNodesInIdOrder`,
  `recalcNodeIdsPostfixOrderAndAddInOrder` (never referenced past its
  own inline definition).
- `TreeNode::getDegree` (used only by the dead `isBinary`/
  `reRootAtHighDegreeNode`/`addNodesWithDegree`/`collapse`),
  `getNumChildren` (zero uses anywhere, even internally),
  `isDescendantOf` (used only by the dead `moveNode`).
- `Tree::drawTreeIDs`/`drawSubtreeIDs` - separate from the structural-
  edit cluster above, but also zero callers anywhere.

A fully commented-out ~127-line block at `Tree_impl.hpp:1495-1621`
(4 complete methods: `assignNodeIdsInUniqueWayFromLeafs`,
`getEulerianWalk`, `getHashCodeFromEulearianWalk`,
`eulerianWalksEquals`) is separately, unambiguously dead - inert
source text, not even compiled.

**Not yet independently re-verified per-method** (the survey above was
reconnaissance to scope this plan, not the same standard of diligence
this engagement applies before actually deleting something) - re-grep
each method individually before removing it, the same way every prior
dead-code phase this engagement has done. Given the cluster's closed
shape, deleting it is likely one phase, but confirm each method's
zero-caller status stands on its own right before deleting, not just
as an inference from cluster membership.

## Finding 2 - `initSubtreeFromStream` calls `exit(1)` on malformed input

`Tree_impl.hpp:150,162` (approximate - re-check exact lines at
execution time): on malformed Newick input, prints `"WRONG INPUT FOR
TREE"` to `std::cerr` and calls `exit(1)` directly, bypassing the
`THROW_EXCEPTION`/`catch`/`exit(EXIT_FAILURE)` convention
`legacy_error_handling_plan.md` established as standard everywhere
else in this codebase (stack unwinding, RAII cleanup, a single
consistent error-reporting path). This is a real, if narrow-trigger,
correctness gap: a library-level parsing routine terminating the
whole process instead of throwing is exactly the anti-pattern that
plan fixed at every other call site it found. Convert both to
`THROW_EXCEPTION`.

## Finding 3 - structural loose ends worth a look, not yet confirmed as bugs

- `TreeNode`'s template copy constructor (`Tree_impl.hpp`, ~35 lines)
  duplicates the non-template copy constructor (~35 lines) almost
  verbatim. Worth understanding why before consolidating - might be a
  genuine necessity of the cross-template-instantiation copy use case
  (`Tree`'s own analogous templated-copy constructor exists for the
  same reason), not necessarily removable.
- `NULL` vs `nullptr` used inconsistently throughout both files (~75
  `NULL` sites in `Tree_impl.hpp` alone, inconsistent even between
  near-twin functions). Purely mechanical, low risk, but real line
  count - worth its own pass rather than folding into every other
  finding's diff.
- Three `(void *)`/C-style casts in assertion/exception-comparison code
  (`Tree_impl.hpp:752,754,788` - line numbers will have shifted after
  the bugfix commit, re-check). Low value to modernize (they're inside
  now-`NDEBUG`-gated `assertTreeStructure`/`THROW_EXCEPTION` checks,
  not hot-path code) - candidate to fold into the `NULL`/`nullptr`
  pass rather than its own phase.
- Three `"(pending const functions)"` comments (`Tree.hpp:264,271,278`)
  flagging unfinished const-correctness on `addInternalNodes`/
  `addNodesWithDegree`/`addNodesOnPath*` - worth a look, but
  `addNodesWithDegree` is only reachable from the dead cluster (Finding
  1), so may become moot once that's deleted.

## Phases (draft - refine once this plan is actually picked up)

1. **Delete Finding 1's dead cluster** (the commented-out block is a
   zero-risk first step within this phase; the ~14-method cluster is
   the real work - re-verify each method individually, delete in one
   pass since the cluster is mutually self-contained, full
   rebuild+test+`RunExamples.sh` after).
2. **Fix Finding 2** (`initSubtreeFromStream`'s `exit(1)` →
   `THROW_EXCEPTION`) - small, independent, verify with a malformed-
   Newick-input manual test in addition to the standard suite (existing
   fixtures likely don't exercise the error path).
3. **Finding 3's loose ends** - `NULL`→`nullptr` mechanical sweep (and
   fold the 3 C-style casts into it), investigate the `TreeNode`
   template-copy-ctor duplication and either consolidate or document
   why it can't be, resolve or drop the 3 pending-const-function
   comments depending on what Phase 1 left reachable.

Each phase gets the standard verification: full rebuild (Clang + real
GCC 14, since these are among the most widely-included headers in the
project) + `ctest` + `RunExamples.sh` byte-identical + `RunCliChecks.sh`,
plus a `-DWITH_LIBXML=OFF` build.

## Explicitly out of scope

- Anything already fixed in commit `f55092f` (`TREE_ASSERT` NDEBUG
  gating, `operator=` self-assignment guard, `makeCanonical`'s
  shadowing bug) - done, not part of this plan.
- A deeper redesign of the intrusive parent/child/sibling linked-list
  structure itself (raw owning pointers, manual `new`/`delete`
  throughout) - consistent with how `Object.hpp`'s removal treated the
  base data structure as out of scope; this plan is about dead code and
  consistency, not an architectural rewrite.
- `fastprot_mpi` - consistent with every other phase in this
  engagement, not verified (no MPI available).
