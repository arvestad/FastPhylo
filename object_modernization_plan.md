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

`Tree::printOn`'s own `os << root` (a `TreeNode*`, relying on
`Object`'s pointer `operator<<` kept alive by the previous commit for
`NJMatrix`'s sake) was left untouched - still compiles unchanged via
that pointer overload for now. Revisit in Phase 4: once `: public
Object` is removed from `TreeNode`, this becomes `os << *root` using
the new reference-form `operator<<` instead.

Verified: full rebuild (Clang + real GCC 14, per this engagement's
standing rule for widely-included headers) + `ctest` 4/4 +
`RunExamples.sh` 20/20 byte-identical + `RunCliChecks.sh`, plus a
`-DWITH_LIBXML=OFF` build - all green, confirming the dispatch-mechanism
change produced byte-identical output for every fixture.

### Phase 2 - de-virtualize `objInitFromStream`; delete `DistanceRow`'s broken constructor

Drop `virtual`/`override` from `objInitFromStream()` on `Sequence`,
`DistanceRow` (well - delete it there, see below), `Tree`,
`DistanceMatrix`, `BitVector`; it becomes a plain method on each,
called exactly the same way it already is (directly, by name, never
polymorphically - Finding 1). Delete `DistanceRow`'s broken
`istream&`-taking constructor and its dead `objInitFromStream()`
declaration/definition at the same time (Finding 3) - it never worked
and has no callers.

This phase is also what removes the `NJMatrix`/`Object*`-pointer-
operator requirement from the previous commit: once `objInitFromStream`/
`printOn` aren't `virtual`, `DistanceMatrix<TreeNode<...>*, ...>`
(`NJMatrix`) no longer eagerly instantiates them for a vtable that no
longer exists - they become ordinary template members, lazily
instantiated only if actually called (which, per the prior commit's
finding, they never are for `NJMatrix`). Expect `Data_init<T*>`/
`Data_printOn<T*>`'s generic `in >> d`/`os << i` to simply never
instantiate for that specialization once this lands - if it does
still get instantiated for some other reason, `TreeNode` will already
have its own `operator<<`/direct-call `objInitFromStream` from this
plan's earlier work, so a pointer-forwarding overload can be added
narrowly rather than resurrecting `Object*`'s generic one.

### Phase 3 - `Tree` equality/hashing; delete the 3 dead `equals()`/`hashCode()` overrides elsewhere

Replace `Tree::equals(const Object*)`/`hashCode()` with a same-typed
`bool operator==(const SequenceTree&, const SequenceTree&)` (or a
`Tree<Data,...>` member, matching how narrowly the real caller -
`tree2int_map` - actually needs it; check whether `Tree<>`'s other
instantiations need this too or whether it can live directly on
`SequenceTree`) and a `std::hash<SequenceTree>` specialization,
replacing `objeq`/`objhash` at `tree2int_map`'s declaration
(`SequenceTree.hpp:110`). Delete `Sequence::equals()`/`hashCode()`,
`BitVector::equals()`/`hashCode()`, `DNA_b128_String::equals()`
outright (Finding 2 - confirmed dead, not being replaced by anything).
Verify `fnj`'s bootstrap-count output is unchanged (`RunExamples.sh`'s
`-b`/bootstrap-flag examples, byte-identical) - this is the one phase
touching actual runtime behavior of a live feature, not just dispatch
mechanism.

### Phase 4 - delete `Object.hpp`/`.cpp`

Once all 7 classes have lost `: public Object`, delete
`Object.hpp`/`.cpp` entirely, including the `operator<<(ostream&,
const Object*)`/`operator>>(istream&, Object*)` pair kept alive by the
previous commit specifically for `NJMatrix`'s sake (should be provably
unnecessary by this point per Phase 2) and the `objeq`/`objhash`
functor structs (superseded by Phase 3's `std::hash`
specialization). Remove the now-dead `#include
"fastphylo/core/Object.hpp"` from all 7 former derived classes' headers
plus any file relying on it transitively for `objhash`/`objeq`
(`SequenceTree.hpp` - already includes `Object.hpp` directly per a
fix earlier this engagement, so update in place rather than remove
outright). Full rebuild (Clang + GCC 14 + `-DWITH_LIBXML=OFF`) +
`ctest` + `RunExamples.sh` + `RunCliChecks.sh` one final time.

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
