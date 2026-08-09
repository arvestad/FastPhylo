# Plan: `Object.hpp`/`.cpp` modernization

## Goal

`Object` is a hand-rolled, Java-style common base (`virtual
toString()`/`printOn()`/`hashCode()`/`equals()`/`operator==`/
`operator!=`/`objInitFromStream()`) inherited by 7 classes across
`core/` and `dna/`. `object_dead_code` cleanup (previous two commits,
2026-08-09) already removed the operators with zero callers
(`objeq_ptr`/`objhash_ptr`/`objeq_ref`/`objhash_ref`, the 4
`operator+` overloads, `operator>>(istream&, Object&)`). What's left
is genuinely reachable code, not dead weight - this plan replaces the
virtual-dispatch design itself with idiomatic modern C++ (free
`operator<<` per class found via ADL, same-typed `operator==`,
`std::hash<T>` specializations), which removes a real latent
undefined-behavior risk along the way (Finding 4) and lets
`Object.hpp`/`.cpp` be deleted outright once nothing derives from it.

Classes deriving from `Object`: `Sequence`, `DistanceRow`, `Tree`,
`TreeNode`, `DistanceMatrix`, `BitVector`, `DNA_b128_String`.

## Finding 1 - what's genuinely load-bearing today, and what isn't

Confirmed via grep for direct/indirect call sites of each virtual
member, not assumed:

- **Printing** (`printOn`/`toString`/`operator<<(ostream&, const
  Object&)`): genuinely live. None of `Sequence`, `Tree`, `TreeNode`,
  `DistanceMatrix`, `BitVector`, `DNA_b128_String` define their own
  `operator<<` - every `os << someSequence`-style call site anywhere
  in the tree resolves through the inherited base via an implicit
  upcast. `DistanceRow` is the exception: it doesn't override
  `printOn` at all, and grep found zero `<<` call sites on any
  `DistanceRow`/`StrDblRow`/`StrFloRow` instance anywhere - it's never
  printed, so it needs no `operator<<` replacement in Phase 1.
- **`equals()`/`hashCode()`/`operator==`/`operator!=`**: only one live
  external path into any of this exists in the whole tree:
  `SequenceTree.hpp`'s `tree2int_map` (`unordered_map<SequenceTree,
  int, objhash, objeq>`), used by `fnj` to count recurring tree
  topologies across bootstrap replicates (`main.cpp`'s
  `buildTrees()`/`tree2count.find(tree)`). That map's key type is
  fixed to `SequenceTree`, so in practice only `Tree<Data,...>`'s
  `equals()`/`hashCode()` override (`Tree_impl.hpp:854`/`876`) is ever
  actually invoked. See Finding 2 for what that means for the other 3
  classes that also override these.
- **`objInitFromStream()`**: overridden by 6 classes (`Sequence`,
  `DistanceRow`, `Tree`, `DistanceMatrix`, `BitVector`, plus `Object`
  itself's no-op default), but never called polymorphically anywhere.
  Every real call site (`DistanceMatrix_impl.hpp:61`,
  `DistanceRow_impl.hpp:40`, `Sequence.cpp:75`'s own body, etc.) calls
  it directly by name on the concrete type from within that class's
  own code. It doesn't need to be `virtual` at all.

## Finding 2 - `Sequence`/`BitVector`/`DNA_b128_String`'s `equals()` (and 2 of the 3's `hashCode()`) overrides are dead code

Grepped every `.equals(`/`.hashCode(`/`operator==`/`operator!=`
call site in the tree with a static type of any `Object`-derived
class. Found exactly one: `objeq`'s own body (`Object.hpp:61`, `o1.equals(o2)`),
which is the `tree2int_map` path from Finding 1. Nothing anywhere
calls `.equals()`/`.hashCode()`/`==`/`!=` directly on a whole
`Sequence`, `BitVector`, or `DNA_b128_String` instance, polymorphically
or otherwise.

That makes `Sequence::equals()`+`hashCode()`
(`Sequence.hpp:52`/`54`, `Sequence.cpp:121`), `BitVector::equals()`+
`hashCode()` (`BitVector.hpp:75`/`76`, `BitVector.cpp:119`), and
`DNA_b128_String::equals()` (`DNA_b128_String.hpp:128`,
`DNA_b128_String.cpp:541`; it doesn't override `hashCode()`) plain
dead code - unlike the `NJMatrix`/`Object*` pointer-operator case from
the previous commit, these are non-template classes, so there's no
vtable-eager-instantiation reason they need to exist. Delete outright
in Phase 3 alongside the `Tree` equality redesign, rather than
"modernizing" logic nothing calls.

## Finding 3 (bonus) - `DistanceRow(std::istream&)` is a broken, self-admittedly non-functional constructor

```cpp
// Doesn't function here :)
DM_TEMPLATE
DISTANCEROW::DistanceRow(std::istream &in)
{
    columns = -1;
    objInitFromStream(in);
}
```

(`DistanceRow_impl.hpp:34-39`, comment included verbatim - the
original author already knew.) Since `DistanceRow` never overrides
`objInitFromStream()`, this calls `Object`'s no-op default
(`{ return is; }`) - nothing is read, and `columns` is left at
`(size_t)-1`. Grep confirms zero callers of this constructor (or of
`DistanceRow`/`StrDblRow`/`StrFloRow` construction from an `istream`)
anywhere in the tree. Not an active bug - just delete it in Phase 2
when `objInitFromStream` stops being an overridable virtual, rather
than leave dead, broken code sitting next to the fix.

## Finding 4 - the actual UB risk being removed

```cpp
TREE_TEMPLATE bool TREE::equals(const Object *o) const
{
    if (getNumNodes() != ((const TREE *)o)->getNumNodes())
```

A raw C-style cast from `Object*` to the concrete template
instantiation `TREE*`, with no `dynamic_cast`/type check. Safe today
only because `tree2int_map`'s key type is pinned to exactly
`SequenceTree`, so both operands `objeq`/`objhash` ever pass through
here really are that type. Nothing in the type system enforces that,
though - `objeq`/`objhash`/`Tree::equals` are declared generically
enough to accept any `Object`, and if either were ever reused for a
differently-typed hashmap, or two different `Tree<>` instantiations
were ever compared through the shared `Object*` interface, the cast
would silently reinterpret one type's memory as another's layout -
undefined behavior per the standard (no guaranteed crash, no
guaranteed error, just unconstrained). A same-typed
`operator==(const Tree&, const Tree&)` (Phase 3) makes this a
compile error instead of a latent runtime hazard - the type mismatch
becomes structurally impossible rather than just "not currently hit."

## Phases

### Phase 1 - printing: per-class free `operator<<`

**Done.** Moved each class's `printOn(ostream&) const` body into a free
`operator<<(std::ostream&, const T&)` found via ADL, for `Sequence`,
`Tree`, `TreeNode`, `DistanceMatrix`, `BitVector`, `DNA_b128_String`
(6 classes - `DistanceRow` needed nothing, confirmed via Finding 1: it
never overrides `printOn` and grep found zero `<<` call sites on it
anywhere). `printOn()` itself stays as a one-line forwarder
(`return os << *this;`) so the class still satisfies `Object`'s virtual
interface unchanged until Phase 4 removes `: public Object` - overload
resolution already prefers the new exact-match free function over the
inherited, conversion-requiring `Object::operator<<` at every real call
site, so `printOn()` becomes unreachable dead code immediately, just
not yet deleted.

`toString()` turned out to have **zero callers anywhere in the tree**,
on any of the 7 `Object`-derived classes (grepped `.toString(`/
`->toString(` tree-wide) - not ported to any class, left as dead
`Object`-inherited code to be removed with the rest of `Object.hpp` in
Phase 4.

Two of the six needed more than a pure cut-and-paste:
- `Tree`/`TreeNode`'s free `operator<<` needed `friend` access inside
  `Tree` (both overloads read `Tree`'s private `dataPrintOn` functor,
  directly for `Tree` and via `TreeNode::getTree()` for `TreeNode`) -
  standard C++ idiom for a free-function `operator<<` needing private
  state, not a design compromise. `TreeNode::printOn`'s internal
  recursive `os << child` (pointer) calls became `os << *child`
  (reference) to match the new signature - safe, since the surrounding
  loop already guarantees `child != nullptr` before each one.
- `DistanceMatrix`'s free `operator<<` needed the same `friend`
  treatment (private `size`/`identifiers`/`D`/`idPrintOn`/
  `distPrintOn`).
- `Sequence`, `BitVector`, `DNA_b128_String` needed no `friend` -
  their `printOn` bodies only ever touched already-public accessors.

**Correction (caught while writing Phase 4, not accurate when first
written):** the paragraph above originally claimed `Tree::printOn`'s
`os << root` was "left untouched," relying on `Object`'s pointer
`operator<<`. That was wrong - the actual code already used the
dereferenced form (`os << *t.root`, `os << *child`) throughout, here
in Phase 1. There was nothing left for Phase 4 to revisit on this
point.

Verified: full rebuild (Clang + real GCC 14, per this engagement's
standing rule for widely-included headers) + `ctest` 4/4 +
`RunExamples.sh` 20/20 byte-identical + `RunCliChecks.sh`, plus a
`-DWITH_LIBXML=OFF` build - all green, confirming the dispatch-mechanism
change produced byte-identical output for every fixture.

### Phase 2 - delete `DistanceRow`'s broken constructor

**Done, smaller than originally scoped.** Deleted `DistanceRow`'s
broken `istream&`-taking constructor and the "Doesn't function here :)"
comment above it (Finding 3) - it never worked (called `Object`'s
no-op default `objInitFromStream`, since `DistanceRow` never overrides
it) and grep confirmed zero callers anywhere.

**Correction to this phase's original text above** (caught while
executing, not before): dropping `virtual`/`override` from
`objInitFromStream()` on `Sequence`/`Tree`/`DistanceMatrix`/
`BitVector` would have been a no-op, not a real de-virtualization.
C++ propagates virtualness from the base automatically whenever the
signature matches and the class still derives from `Object` -
removing the derived declaration's own `override` keyword doesn't
remove the vtable entry; it only removes the compiler's override-
mismatch safety check, for zero behavioral benefit. The vtable entry
only actually disappears once `: public Object` itself is removed,
which doesn't happen until Phase 4. Moved that part there instead of
doing something cosmetic-only now.

Same correction applies to this phase's `NJMatrix`/`Object*`-pointer-
operator claim: that requirement isn't removed by touching
`objInitFromStream`'s keywords at all - it goes away in Phase 4, when
`DistanceMatrix`/`TreeNode` actually stop deriving from `Object` and
`Data_init<T*>`/`Data_printOn<T*>` have no more virtual vtable to
eagerly instantiate for. Re-verify this claim concretely there rather
than assuming it from here.

Verified: full rebuild (Clang + GCC 14) + `ctest` 4/4 + `RunExamples.sh`
byte-identical + `RunCliChecks.sh`.

### Phase 3 - `Tree` equality/hashing; delete the 3 dead `equals()`/`hashCode()` overrides elsewhere

**Done.** Replaced `Tree::equals(const Object*)`/`hashCode()` with a
free, same-typed `operator==`/`operator!=(const Tree<Data,...>&, const
Tree<Data,...>&)` template (kept generic at the `Tree<>` level, same
scope the old virtual overrides had - only `SequenceTree` has a live
caller today, but nothing else needed narrowing it) and a
`std::hash<Tree<Data,DataInitializer,DataPrintOn>>` partial
specialization. Updated `tree2int_map` (`SequenceTree.hpp`) to
`std::unordered_map<SequenceTree, int>`, relying on the new operators
via their default template arguments instead of the `objeq`/`objhash`
functors. Deleted `Sequence::equals()`/`hashCode()` (+ the now-orphaned
`stringhasher` static member and its `hashstr`/`stl_utils.hpp`
dependency), `BitVector::equals()`/`hashCode()`, and
`DNA_b128_String::equals()` (+ its only caller, a backwards
`ambiguity_nucleotide_at_position` `operator==` that returned `true`
for *different* values - dead code, not a live bug, deleted rather
than fixed) - all confirmed dead per Finding 2.

**One real subtlety, not anticipated in this plan's original text:**
`SequenceTree` is a distinct class that *derives from*
`Tree<Sequence_double, ...>`, not an alias for it. Template partial
specialization matching is structural/exact - `std::hash<Tree<Data,
...>>` does **not** cover `SequenceTree` on its own, unlike
`operator==`, where ordinary template argument deduction succeeds
through a derived-to-unique-base relationship. Needed one additional
explicit (non-partial) specialization, `std::hash<SequenceTree>` in
`SequenceTree.hpp`, delegating to the `Tree<>` one via the same
implicit upcast `operator==` already relies on. Caught by the Clang
build failing outright (`hash<SequenceTree>`'s implicitly-deleted
default constructor) - not a silent issue.

**A second real bug caught only by the mandatory real-GCC-14 rebuild**
(not Clang): removing `stl_utils.hpp` from `Sequence.hpp` (no longer
needed once `stringhasher`/`hashstr` were gone) also dropped a
transitive `<vector>` include that `Sequence.hpp` itself needs for its
own `std::vector<Sequence>` parameters - invisible on Clang/libc++
(something else in the chain still pulled it in) but a hard compile
error on GCC/libstdc++. Fixed by adding `#include <vector>` directly
to `Sequence.hpp`, the file that actually uses it - same category of
bug, and same fix, as the `<iomanip>` case earlier this engagement.

Verified: full rebuild (Clang + real GCC 14) + `ctest` 4/4 +
`RunExamples.sh` 20/20 byte-identical (including example 19's `-b 3`
bootstrap pipeline, the one real live consumer of this machinery) +
`RunCliChecks.sh`, plus a manual `fnj --print-counts` run confirming
bootstrap replicate counting is still correct, plus a
`-DWITH_LIBXML=OFF` build.

### Phase 4 - delete `Object.hpp`/`.cpp`

**Done, and larger than planned: `objInitFromStream` turned out to be
fully dead, not merely non-polymorphic.** Before touching anything,
re-verified Phase 1's Finding 1 claim ("`objInitFromStream` is called
directly, by name, never polymorphically") rather than trusting it.
It was wrong for `Sequence`: zero callers anywhere, direct or
otherwise - `Sequence::readSequences` has its own separate parsing
logic, never calling it. Same story for `Tree`: its own
`objInitFromStream` body doesn't even call itself recursively for
anything real, and grepping for constructions of any `Tree<>`
instantiation via `istream` (`Tree<...>(someIstream)`,
`SequenceTree(someIstream)`) outside the class's own declaration files
found zero call sites (`SequenceBasedNJ.cpp`'s `b128Tree` typedef is
the only other `Tree<>` instantiation anywhere, and that file has been
unbuilt for multiple phases now). Same for `DistanceMatrix`'s
`objInitFromStream` and its `istream&`-taking constructor - zero
callers, confirmed by grep and then by the compiler itself (the
rebuild succeeded through `DistanceMatrix.cpp` with the constructor
gone, without complaint). `BitVector`'s `objInitFromStream` override
was an inline no-op with no `istream&` constructor to call it at all.

**Deleted outright, not just de-virtualized:** `objInitFromStream` on
`Sequence`/`Tree`/`DistanceMatrix`/`BitVector` (declarations and
definitions), `DistanceMatrix`'s now-pointless `istream&` constructor,
and `DistanceMatrix`'s now-orphaned private `idInit`/`distInit`
functor members (used nowhere else - `idPrintOn`/`distPrintOn` stay,
genuinely used by Phase 1's `operator<<`). Also deleted every
`printOn()` forwarder across all 6 classes Phase 1 touched (Sequence,
Tree, TreeNode, DistanceMatrix, BitVector, DNA_b128_String) - Phase 1
already established nothing calls them; Phase 4 is what finally
removes the dead weight rather than leaving a forwarder with no
purpose. `~Tree()`/`~TreeNode()` changed from `override` to plain
`virtual` (not removed - `SequenceTree` is a real, live subclass of
`Tree<Sequence_double,...>`, so polymorphic destruction through a
`Tree<>*`/`TreeNode*` base pointer must stay safe regardless of
`Object`; none of the other 5 classes have subclasses, so their
implicit destructors needed no change at all).

**A real mistake, caught by the build, not by review:** initially
deleted `Tree(std::istream&)`/`SequenceTree(std::istream&)` too,
believing them dead by the same reasoning as `objInitFromStream` (which
that constructor doesn't even call - it goes through the separate,
very much live `initSubtreeFromStream`, also used by the Newick-string
constructor). They aren't dead:
`tests/SequenceTree_PhylipReader_test.cpp:44` builds a tree from an
`istringstream` via exactly this constructor. Missed on the first grep
pass because the check looked for the literal call pattern
(`SequenceTree(someIstream)`) but didn't verify what type `newick`
actually was at that call site - it read as a Newick *string* by
name, but is in fact declared as `istringstream newick(...)` two lines
above. Restored both constructors (verified independent of
`objInitFromStream`, so no conflict with deleting that); the failing
`ctest` run caught it immediately, before any commit.

**The `NJMatrix`/`Object*`-pointer-operator concern resolved itself
exactly as predicted**, and was concretely re-verified rather than
assumed: once `TreeNode`/`DistanceMatrix` no longer derive from
`Object`, `NJMatrix` (`DistanceMatrix<TreeNode<...>*, ...>`) has no
vtable left to eagerly instantiate `objInitFromStream`/`printOn`
for, so `Data_init<TreeNode*>`/`Data_printOn<TreeNode*>`'s generic
`in >> d`/`os << i` never got instantiated at all - the build went
straight through with no pointer-operator overload needed anywhere.
(It also turned out Phase 1's own `Tree`/`TreeNode` `operator<<`
already used the dereferenced form throughout, not the pointer form a
misremembered note in this plan's Phase 1 writeup claimed was "left
untouched" - that note was simply inaccurate, corrected here.)

Deleted `Object.hpp`/`.cpp` entirely, including the
`operator<<(ostream&, const Object*)`/`operator>>(istream&, Object*)`
pair from the dead-code-cleanup commit (confirmed unnecessary per
above) and the `objeq`/`objhash` functor structs (superseded by Phase
3's `std::hash` specialization). Removed the `#include
"fastphylo/core/Object.hpp"` from all 7 former derived classes'
headers and the `Object.cpp` entry from `src/c++/CMakeLists.txt`.
Grepped the whole tree for any remaining reference to `Object`/
`objhash`/`objeq` afterward - none outside historical prose comments
and the plan doc itself.

Verified: a full clean rebuild from scratch (removed `build/` entirely
first, not just incremental) + Clang + real GCC 14 + `ctest` 4/4 +
`RunExamples.sh` 20/20 byte-identical + `RunCliChecks.sh`, plus a
`-DWITH_LIBXML=OFF` build - all green.

## Explicitly out of scope

- `fastprot_mpi` - consistent with every other phase in this
  engagement, not verified (no MPI available); changes applied there
  by inspection/analogy only.
- Any behavior change to what gets printed/compared/hashed - this is
  strictly a dispatch-mechanism change (virtual → free
  function/same-type), not a reformatting or feature change. Phase 3
  is the only phase touching a live feature's code path
  (`tree2int_map`), and even there the goal is identical output, just
  produced through a type-safe comparison instead of an unchecked
  cast.
- `Tree.hpp`/`Tree_impl.hpp`'s broader audit - out of scope per
  `modern_cpp_core_cleanup_plan.md` already; this plan only touches
  the specific `Object`-derived interface (`printOn`/`equals`/
  `hashCode`/`objInitFromStream`), not the rest of those files' 2000+
  lines.
