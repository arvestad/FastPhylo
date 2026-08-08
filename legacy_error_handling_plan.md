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

## Phases

### Phase 0 - remove the dead, unbuilt legacy apps (Finding 4)

Confirm none of `sequence_nj.cpp`/`buildtree.cpp`/
`iterative_tree_merger/`/`tests/DistanceMatrix_test.cpp` is referenced
by any build system, doc, or script anywhere (grep already suggests
not); delete them. Independent of every other phase - do this first,
it's the same kind of call already made for `Simulator.hpp/cpp` and
removes dead weight before it can complicate later greps.

### Phase 1 - fix the exit-code bug (Finding 0)

Add `exit(EXIT_FAILURE);` to both catch branches in all four apps'
`main()`. Trivial, high-value, zero design risk - do this regardless
of how the rest of the plan goes, and independently of Phase 0.

### Phase 2 - replace the macros with real C++ (Finding 1)

Per app: replace `TRY_EXCEPTION(); ... CATCH_EXCEPTION() catch(...){...}`
with a plain, visible `try { ... } catch (const Exception& e) { ... }
catch (const std::exception& e) { ... }` (now that Finding 3 fixes
`what()`, the second catch becomes genuinely useful instead of a
`std::terminate`-avoidance stopgap). Delete `TRY_EXCEPTION`,
`CATCH_EXCEPTION`, `CATCH_RETHROW` (already confirmed dead) from
`Exception.hpp`. Decide whether `THROW_EXCEPTION` stays as-is (it
still pulls its weight for `operator<<`-style messages) or becomes a
small variadic helper - lower priority, can go either way without
blocking the rest.

### Phase 3 - fix `Exception`'s relationship to `std::exception` (Finding 2)

Add `const char* what() const noexcept override` returning `message`
(or a cached formatted string, if `printOn`'s output should be
preserved verbatim). Decide whether `Exception` still needs to inherit
`Object` - grep how much of `Object`'s interface (`hashCode`/`equals`/
`toString`) `Exception` actually uses beyond `printOn`; if it's just
`printOn`, that could be `Exception`'s own method instead, dropping
the `Object` dependency entirely and leaving a small, standard-shaped
exception type.

### Phase 4 - modernize `file_utils` (Finding 3)

Delete the 7 dead functions. Replace `file_exists` with
`std::filesystem::exists` at its call site(s) and delete the wrapper.
Convert the remaining live functions' primary signatures to
`std::string`/`std::string_view`/`std::filesystem::path` as
appropriate, collapsing each function to one real implementation
rather than a `char*` original plus a `.c_str()`-forwarding
`std::string` overload.

## Explicitly out of scope

- `fastprot_mpi` - consistent with every other phase in this
  engagement, not verified (no MPI available); apply the same changes
  by inspection/analogy only if asked.
- Any change to `Object.hpp` itself beyond what Phase 3 needs from
  `Exception`'s side - `Object` is used well beyond `Exception`
  (`Sequence`, `DistanceMatrix`, etc.); a full `Object` redesign is a
  separate, much bigger decision not raised here.
