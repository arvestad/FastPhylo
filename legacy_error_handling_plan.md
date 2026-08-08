# Plan: legacy error-handling and file-I/O cleanup

## Goal

`Exception.hpp`/`file_utils.hpp` are pre-C++11, Java-flavored
infrastructure from the original 2006 codebase that no prior pass
fully addressed - `lint_plan.md`'s Phase 1 patched one symptom (the
`bugprone-exception-escape` gap, via an outer `catch(...)`) without
touching the underlying macros, the `Exception` class design, or
`file_utils`'s dead functions. This plan finishes that job: modernize
*and* simplify, per direct instruction - not just make it compile
cleanly, actually reduce how much of this a reader has to hold in
their head.

Found while investigating (2026-08-08), verified via grep/read, not
assumed:

## Finding 0 - a live correctness bug: apps exit 0 on error

All four apps' `main()` (`fastdist`, `fnj`, `fastprot`, and by the same
pattern presumably `fastprot_mpi`) end with:

```cpp
try {
  ...
}
catch(...){
  throw;
}
CATCH_EXCEPTION()          // }catch(Exception exc){ std::cerr << exc << std::endl; }
catch(...){
  std::cerr << "Unknown (non-Exception) error" << std::endl;
}
return 0;
```

Neither catch branch calls `exit(EXIT_FAILURE)` - both fall through to
the unconditional `return 0;`. **An app that hits a real error prints
a message to stderr and still exits with status 0 (success).** Any
script or pipeline checking the exit code is silently lied to. This is
the highest-priority, lowest-risk item here: a one-line fix
(`exit(EXIT_FAILURE);` at the end of each catch branch) per app, no
design decisions needed. Should happen first, independent of
everything else below.

## Finding 1 - `TRY_EXCEPTION`/`CATCH_EXCEPTION`/`CATCH_RETHROW` macros

`include/fastphylo/core/Exception.hpp`:

```cpp
#define THROW_EXCEPTION(MES)  {std::ostringstream out; out << MES; throw Exception(__FILE__,__FUNCTION__,__LINE__,out.str());}
#define TRY_EXCEPTION() try{
#define CATCH_EXCEPTION() }catch(Exception exc){ std::cerr << exc <<std::endl;}
#define CATCH_RETHROW() }catch(Exception exc1){ Exception exc2(__FILE__,__FUNCTION__,__LINE__,""); exc2.addToStackTrace(exc1); throw exc2;}
```

`TRY_EXCEPTION()`/`CATCH_EXCEPTION()` each expand to *half* a
`try`/`catch` - they only produce valid C++ because each call site
happens to balance the braces textually across two macro invocations
with arbitrary code in between. No editor, formatter, or reader can
see the real control flow at the macro call site. `CATCH_EXCEPTION()`
also catches `Exception exc` **by value** - an unnecessary copy on
every exception, and exactly what `clang-tidy`'s
`misc-throw-by-value-catch-by-reference` flags.

`THROW_EXCEPTION` itself (a normal object-like function macro building
an `ostringstream` and throwing) is more defensible - it exists so
call sites can write `THROW_EXCEPTION("bad value: " << x)` with
`operator<<`-style message building, which a plain function can't do
as tersely. Worth keeping *if* still needed after Finding 2, or
replacing with a small variadic/stream helper - see Phase 2.

`CATCH_RETHROW()` - grep confirms zero call sites anywhere in the
tree. Dead, delete outright regardless of what else happens here.

## Finding 2 - the `Exception` class double-inherits, and it's broken

```cpp
class Exception : public Object, public std::exception
```

`Object` (`include/fastphylo/core/Object.hpp`) is a hand-rolled,
Java-style base (`toString()`/`hashCode()`/`equals()`/`operator==`)
used across the old core classes. `Exception` also inherits
`std::exception` - presumably so it's catchable via the standard
`catch (const std::exception&)` idiom - but **never overrides
`what()`**. Anyone who did catch it as `const std::exception&` and
called `.what()` would get the useless default `"std::exception"`
string; the actual message only surfaces through `Exception`'s own
`printOn()`/`operator<<`. Two unrelated error-reporting conventions
bolted together, neither fully realized.

## Finding 3 - `file_utils.hpp`: half the functions are dead, all of it is `char*`-based

**Superseded by Phase 4's re-check below** - this table was correct at
the time (2026-08-08, before Phase 0 landed) but Phase 0 deleted the
only remaining callers of two of the "live" entries here, and two more
turned out to be mis-assessed (an overload-name collision with
`stl_utils.hpp`). Kept for the historical record of the original
investigation; see Phase 4 for the final, correct live/dead count (3
of 14 survive, not 7).

Grepped every declared function for call sites outside
`file_utils.cpp`/`.hpp` themselves. 7 of 14 have **zero**:

| Dead (delete) | Live (keep, modernize) |
|---|---|
| `open_write_file_interactive` | `file_exists` |
| `open_write_stream_interactive` | `open_write_file` |
| `open_read_file` | `open_write_stream` (both overloads) |
| `open_read_file_interactive` | `open_read_stream` (both overloads) |
| `open_read_stream_interactive` | `open_write_binary` |
| `skipUntil` | `skipWhiteSpace` (both overloads) |
| `appendToken` | `appendUntil` |

The "interactive" family (prompt-the-user-and-retry on missing/
existing files) is a natural set to lose together - it's a CLI-prompt
UX pattern from a much older, non-scripted usage model, unused by any
of the four current apps.

Direct answer to "is there an actual need for `file_exists(char*)`":
no. This project already requires C++17 (this whole branch is called
`modernize-cpp17`), and nothing in the tree uses `<filesystem>` yet.
`file_exists` should become a direct call to `std::filesystem::exists`
at its 1 call site, or every remaining caller can just call
`std::filesystem::exists` themselves and the function can be deleted
entirely - it's a 4-line wrapper adding no behavior of its own beyond
what the standard library now provides natively.

And yes, on `char*` vs `std::string`: every one of these signatures
predates `std::string` being idiomatic (the file header comments date
to 2006). The 4 functions that already have a `std::string` overload
(`open_write_stream`, `open_read_stream`) only got it as a thin
`.c_str()`-forwarding wrapper around the "real" `char*` version -
backwards from where the codebase should end up. Modernizing this
properly means the `std::string`/`std::string_view` signature becomes
the one real implementation, `FILE*`-based functions take
`const std::filesystem::path&` where they're about paths rather than
in-memory buffers, and the `char*` overloads either go away (implicit
`std::string` construction covers callers passing a string literal for
free) or become trivial one-line forwarders if a `const char*` call
site is found that matters.

## Finding 4 (bonus, found while checking `file_utils`' callers) - a batch of dead, unbuilt legacy apps, one containing a real signature-mismatch bug

Grepping `open_read_stream`'s callers surfaced `src/c++/sequence_nj.cpp`,
`src/c++/buildtree.cpp`, `src/c++/iterative_tree_merger/*.cpp`, and
`tests/DistanceMatrix_test.cpp`. **None of these are referenced by any
`CMakeLists.txt` anywhere in the tree** - confirmed via grep, same
situation as the `Simulator.hpp/cpp` code already deleted earlier in
this engagement (see `modernization_plan.md`). They never compile, so
nothing catches bugs in them.

One such bug: `file_utils.cpp`'s `open_read_stream(const string,
ofstream&)` (line ~212) takes an **`ofstream&`**, but
`file_utils.hpp` declares that overload as `open_read_stream(const
std::string, std::ifstream&)` (`ifstream`, not `ofstream`). Header and
definition disagree - this would fail to *link* (mismatched mangled
symbol) the moment anything called the header's declared overload.
Nothing does, because nothing that calls it is built. This isn't
something to fix - it's a symptom that these files should be deleted,
the same call made for `Simulator.hpp/cpp`. Recommend confirming with
Lasse this code has no future use (same as that prior removal) and
deleting `sequence_nj.cpp`, `buildtree.cpp`,
`iterative_tree_merger/`, and `tests/DistanceMatrix_test.cpp` in
Phase 0 of this plan, before touching `file_utils` itself - it shrinks
the surface area for Finding 3's cleanup and removes the only place
the header/definition mismatch could ever have mattered.

## Finding 5 - full audit of `src/fastphylo/core`/`include/fastphylo/core` (2026-08-08)

Direct instruction: scrutinize every file in `core/`, not just
`Exception`/`file_utils`. Systematic grep sweep of all 34 files (~7200
lines) plus the closely-related `src/fastphylo/dna/LeastSquaresFit.*`.

**More dead code, same shape as Finding 4 (confirmed via grep, zero
live callers):**

- `arg_utils.c`/`.h` + `arg_utils_ext.cpp`/`.hpp` - a pre-CLI11
  leftover (`gengetopt_migration_plan.md`), zero callers anywhere.
  `arg_utils_ext.cpp`/`arg_utils.c` were still listed in
  `src/c++/CMakeLists.txt` and linked straight into `libfastphylo.a` -
  dead code shipping in the binary, not just dead source.
- `SequenceTree_MostParsimonious.cpp` (214 lines, also linked into
  `libfastphylo.a`) - its only callers anywhere in the tree were two
  more dead test files (`Tree_test.cpp`, `ParsimonyTest.cpp`).
  Unreachable from any of the 4 built apps.
- 6 more unbuilt test files beyond `DistanceMatrix_test.cpp` (already
  deleted in Phase 0): `AML_star_test.cpp`, `Big_AML_test.cpp`,
  `DNA_b128_String_test.cpp`, `LeastSquaresFit_test.cpp`,
  `ParsimonyTest.cpp`, `Tree_test.cpp` - none registered in
  `CMakeLists.txt`/`ctest`. The `AML_*` ones test code (`Big_AML.cpp`
  etc.) that doesn't exist in the tree anymore at all.
- `src/fastphylo/dna/LeastSquaresFit.cpp`/`.hpp` - one directory over
  from `core/`, same pattern: linked into the library, zero callers
  from any app.
- `log_utils.hpp`'s `ASSERT`/`ASSERT_OP`/`DEBUG` macros - all zero
  call sites. `DEBUG`'s body has a latent syntax bug (a string literal
  directly adjacent to `__LINE__` with no `<<` between them - would
  fail to compile if ever invoked) and `ASSERT_OP`'s `NDEBUG` branch
  has a copy-paste bug (redefines `ASSERT_EQ` instead of `ASSERT_OP`) -
  both dormant only because nothing calls either macro.

All of the above deleted directly (Phase 5), verified via clean
rebuild + `ctest` 4/4 + `RunExamples.sh` 20/20 byte-identical +
`RunCliChecks.sh` all green after each step.

**Checked and turned out fine - no action needed:** `Object`'s
`equals()`/`hashCode()`/`toString()`/`objInitFromStream()` are all
genuinely used elsewhere in the tree - `equals()`/`hashCode()` via the
`objeq`/`objhash` functor structs (used by `SequenceTree`'s hash
containers), `objInitFromStream()` overridden and called across 6
classes (`Sequence`, `DistanceRow`, `DistanceMatrix`,
`FloatDistanceMatrix`, `Tree`, `BitVector`). `Object` is not dead
weight, just old-style (Java-flavored) C++ - only `Exception`'s
specific double-inheritance from both `Object` and `std::exception`
(Finding 2) is actually broken. No case for touching the other 7
`Object`-derived classes here.

**Real architectural finding, deliberately not acted on - too large
for this plan:** there are **two parallel, uncoordinated
error-handling systems** in this codebase. `Exception`/
`THROW_EXCEPTION` throws a C++ exception, unwinds the stack, and (as
of Phase 1) is caught cleanly in each `main()` and exits 1. But
`log_utils.hpp`'s `USER_ERROR`/`PROG_ERROR`/`MEM_CHECK` macros (33+
call sites across the tree, well beyond `core/`) call `exit(1)`/
`std::exit(1)` **directly** from wherever they're invoked - no stack
unwinding, no RAII cleanup, no chance for any `main()`'s catch blocks
to run. Reconciling these into one consistent error-reporting
convention is a genuine design decision (which one wins? does
`USER_ERROR` become `THROW_EXCEPTION` everywhere, changing behavior
at 20+ call sites?) touching far more of the tree than this plan's
scope - flagging for a dedicated future plan, not attempting here.

## Phases

### Phase 0 - remove the dead, unbuilt legacy apps (Finding 4)

**Done** (`deef6af`, 2026-08-08).

Confirm none of `sequence_nj.cpp`/`buildtree.cpp`/
`iterative_tree_merger/`/`tests/DistanceMatrix_test.cpp` is referenced
by any build system, doc, or script anywhere (grep already suggests
not); delete them. Independent of every other phase - do this first,
it's the same kind of call already made for `Simulator.hpp/cpp` and
removes dead weight before it can complicate later greps.

### Phase 1 - fix the exit-code bug (Finding 0)

**Done** (`ba6d676`, 2026-08-08). Fixed at the source rather than per
app: `CATCH_EXCEPTION()` itself now calls `exit(EXIT_FAILURE)` after
printing, covering all four apps (including `fastprot_mpi`, which only
uses that macro and has no separate `catch(...)` fallback) from one
place. The three apps with their own `catch(...)` "Unknown
(non-Exception) error" fallback (`fastdist`/`fnj`/`fastprot`) each got
`exit(EXIT_FAILURE)` added there too.

### Phase 2 + Phase 3 - delete `Exception` entirely, use `std::runtime_error` (Findings 1 + 2)

**Done, superseding the original two-phase plan below.** Original plan
was: fix `Exception`'s `what()` (Phase 3), then replace the macros with
plain `try`/`catch (const Exception&)`/`catch (const std::exception&)`
(Phase 2). The `what()` fix was implemented and confirmed working
(catching as `const std::exception&` now returned the real message
instead of the useless default) - but a direct question ("why isn't
the standard exception good enough?") led to checking whether
`Exception` earns its keep at all, not just whether its bug was fixed.

It doesn't. Verified via grep: nothing anywhere in the tree reads
`Exception`'s `.file`/`.function`/`.line`/`.stackTrace` fields directly,
nothing catches it structurally (by type) except the macros themselves,
and `addToStackTrace`/`CATCH_RETHROW` are entirely dead (Finding 1).
All 57 `THROW_EXCEPTION` call sites across 16 files use identical
generic streaming syntax (`THROW_EXCEPTION("msg" << var)`) - since
it's a macro, none needed to change when its definition changed.

Final shape: `THROW_EXCEPTION(MES)` now builds a `file:line (function):
message` string and throws `std::runtime_error` directly - same
diagnostic info as before (file/line/function), just as message text
instead of struct fields, no custom class. `Exception.hpp`/`.cpp`
deleted outright (Object dependency for this path gone with it).
`TRY_EXCEPTION`/`CATCH_EXCEPTION`/`CATCH_RETHROW` deleted - each app's
`main()` now has a plain, visible `try { ... } catch (const
std::exception& e) { std::cerr << e.what() << ...; exit(EXIT_FAILURE);
} catch (...) { ... }` (the inner no-op `catch(...){throw;}` wrapper
each app had was also removed - it had zero behavioral effect).

One real bug surfaced by this: `stl_utils.hpp` used `Object`/`objhash`/
`objeq` without including `Object.hpp` directly - it relied on
`Exception.hpp` transitively pulling it in. Fixed with a direct
`#include "fastphylo/core/Object.hpp"` rather than restoring the
transitive path.

`fastprot_mpi` updated by inspection only (unverified, no MPI
available) - same mechanical transformation, plus added the
`catch(...)` "Unknown error" fallback it never had (bringing it in
line with the other three apps).

Verified: full rebuild + `ctest` 4/4 + `RunExamples.sh` 20/20
byte-identical + `RunCliChecks.sh` all green, both with the default
build and a separate `-DWITH_LIBXML=OFF` build (to catch any
`Object.hpp`-transitive-include fallout in that code path too). Spot
checked live error output - `fastdist` on a missing file now prints a
single line (`.../PhylipMaInputStream.cpp:22 (PhylipSequenceReader):
File doesn't exist: "..."`) and exits 1, replacing the old 7-line
boxed format - strictly less code, same information.

### Phase 4 - modernize `file_utils` (Finding 3)

**Done, with a larger scope than Finding 3's original table.**
Finding 3 was assessed *before* Phase 0 deleted `sequence_nj.cpp`/
`buildtree.cpp`/`iterative_tree_merger/` - which turned out to be the
only remaining callers of `open_write_stream`/`open_read_stream`.
Re-checking at Phase 4 time (always re-verify current state rather
than trusting an earlier finding once other phases have landed)
found both now fully dead too, plus two more the original table
mis-assessed: `appendUntil` in `file_utils.hpp` (the 3-arg
`istream`-based overload) turned out to share a name with a
completely different 4-arg `appendUntil` in `stl_utils.hpp` - the
original grep matched call sites for the latter, not the former,
which has zero real callers. Likewise `skipWhiteSpace(FILE*)` - only
the `istream&` overload is ever called.

Final live/dead count: of the original 14 declared functions (18
counting overloads separately), only 3 survive:
`open_write_file`, `open_write_binary`, `skipWhiteSpace(std::istream&)`.
Everything else deleted, including `file_exists` (its only 2 callers
were the two already-dead "interactive" functions, so nothing needed
a `std::filesystem::exists` replacement - it just became orphaned and
was deleted directly).

`open_write_file`/`open_write_binary` kept their `const char*`
signature rather than moving to `std::string`/`path` as originally
planned - their only remaining callers (`DataOutputStream`'s
constructor family) pass a nullable `char *filename` where `nullptr`
means "write to stdout", an idiom `std::string`/`string_view` can't
represent without introducing `std::optional` or an overload set -
that's a real design step belonging to `DataOutputStream`, not a
`file_utils.hpp` signature-modernization detail. Not attempted here.

Real bug surfaced while rebuilding: `file_utils.hpp` used to
`#include "fastphylo/core/Exception.hpp"`, and 9 other files across
the tree used `THROW_EXCEPTION` while relying on that entirely
transitively (never including `Exception.hpp` directly themselves).
Once the rewritten `file_utils.hpp` dropped that include (no longer
needed by anything declared in the header itself), all 9 broke. Fixed
correctly - each of the 9 now includes `Exception.hpp` directly -
rather than papering over it by restoring the transitive path.

Verified: full rebuild + `ctest` 4/4 + `RunExamples.sh` 20/20
byte-identical + `RunCliChecks.sh` all green, plus a separate
`-DWITH_LIBXML=OFF` build (several of the 9 broken files were
XML-related).

### Phase 5 - delete the rest of the dead code found in the full `core/` audit (Finding 5)

**Done**, same commit as Finding 5's discovery. `arg_utils`/
`arg_utils_ext` (4 files), `SequenceTree_MostParsimonious.cpp` (+ its
declaration in `SequenceTree.hpp`), the 6 more dead test files,
`LeastSquaresFit.cpp`/`.hpp`, and `log_utils.hpp`'s dead
`ASSERT`/`ASSERT_OP`/`DEBUG` macros - all deleted, `CMakeLists.txt`
updated to match, verified via full rebuild + test suite.

## Explicitly out of scope

- `fastprot_mpi` - consistent with every other phase in this
  engagement, not verified (no MPI available); Phase 2+3's changes were
  applied there by inspection/analogy only (done).
- Any change to `Object.hpp` itself - `Exception` no longer depends on
  it at all (Phase 2+3 deleted `Exception` entirely), and `Object` is
  used well beyond that (`Sequence`, `DistanceMatrix`, etc.); a full
  `Object` redesign is a separate, much bigger decision not raised
  here. Finding 5 confirmed `Object`'s interface is all genuinely used,
  so there's no dead-code case for touching it either.
- Reconciling `Exception`/`THROW_EXCEPTION` with `log_utils.hpp`'s
  `USER_ERROR`/`PROG_ERROR`/`MEM_CHECK` (Finding 5) - a real, separate
  design decision affecting 30+ call sites across the tree, not a
  mechanical cleanup like the rest of this plan.
