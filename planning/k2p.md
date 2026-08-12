K2p

## Resolution (2026-08-12, DONE)

Implemented the recommended fix: swapped `fastdist`'s default for
**both K2P and TN93** (confirmed structurally, not just suspected -
`fillMatrix_TN93` has the identical `no_tstvratio` dispatch shape as
K2P, sharing the same `sequence_translation_model` struct field).

`--tstvratio`/`--pyrtvratio` are now what opts into the fixed-ratio
maximum-likelihood method (using the given value(s)); with neither
given, `fastdist` estimates the ratio per pair (the classical/textbook
formula) by default. `--no-tstvratio` is kept, unchanged in meaning,
for scripts that already pass it - it's a no-op against the new
default now, and errors if combined with `--tstvratio`/`--pyrtvratio`
(previously-contradictory intent, now rejected instead of silently
picking one).

Verified against this doc's own cross-validation methodology: the new
default output on a synthetic dataset matched a hand-computed classical
K2P reference exactly (previously matched only with `-N` explicitly
given); `fastdist -T 2` reproduces the *old* default's fixed-ratio-at-
2.0 numbers exactly, confirming the old path is still reachable, just
opt-in now. `ctest` 7/7 (Clang + real GCC 14), `RunCliChecks.sh` (added
a check for the new `-N`/`-T` conflict rejection), and
`RunExamples.sh` all clean - 3 fixtures (`ex1`/`ex2`/`ex3`, the ones
using bare `fastdist` with no `-D`/`-T`/`-N`) legitimately changed and
were regenerated, since they now correctly reflect the classical
formula instead of the old fixed-ratio-at-2.0 default.

Also found while checking TN93's shared dispatch shape:
`TamuraNei.cpp`'s own `compute_Tamura_Nei_fixratio()` calls a *second*,
separate `secant_search()` (different file from Kimura2parameter.cpp's)
that reads its inputs (`purine_ratio`, `tU`, `tY`, `tV`, `n`,
`gA`/`gC`/`gG`/`gT`/`gU`/`gY`) from file-scope mutable globals instead
of parameters, rather than the bare-float unbounded search
Kimura2parameter.cpp's own `secant_search()` had before its 2026-08-12
fix (unbounded step size that could overflow `exp(2*x*K)` and either
crash to nan or falsely "converge" to a wildly wrong finite value -
see that commit's message for the full story). **Not audited or fixed
here** - it's a distinct robustness question from this doc's own scope
(the wrong-*default* bias, not a numerical-search bug), and the global-
mutable-state shape alone deserves its own careful look before touching
the math. Flagged for a follow-up, not silently left unmentioned.

## Goal

`fastdist -D K2P` (the default distance function - no `-D` needed to
trigger this) computes systematically **inflated** distances relative
to the standard, textbook Kimura (1980) K2P formula, and the inflation
**grows with true divergence**. This was found and confirmed by an
external, downstream consumer of the library
(`dnctree_testing/experiments/alisim_dna_matrix`, a phylogenetic tool
benchmarking project), not from FastPhylo's own test suite, and it
materially affected that project's results before they worked around
it - see "How this was found" below.

The existing `-N`/`--no-tstvratio` flag already fixes this (confirmed
by direct testing), but it isn't the default, its help text doesn't
explain what the default actually does differently, and it's easy to
miss entirely if you're just calling `fastdist -D K2P` expecting "the"
K2P distance.

## How this was found

`dnctree_testing`'s `alisim_dna_matrix` experiment compares several
tree-inference tools on HKY-simulated DNA data, several of which are
fed a `fastdist -D K2P`-computed distance matrix (RapidNJ, BioNJ,
DIPPER) since none of the tools implement true HKY and K2P is the
closest available approximation. One tool in the comparison,
`RapidNJ (standalone)`, computes its **own** K2P distances directly
from the alignment (`rapidnj -a kim -t d`) instead of consuming
`fastdist`'s matrix - and it consistently, meaningfully outperformed
the `fastdist`-fed RapidNJ run on identical data, which is suspicious:
same tool, same tree-building algorithm, the only thing that differs
is the distance *source*.

Direct investigation on one real test alignment (200 taxa,
HKY-simulated, 1500nt) traced this to `fastdist` itself, by comparing
four independent computations of the same specific taxon pairs:

1. The classical closed-form Kimura (1980) K2P formula, hand-computed
   directly from each pair's observed transition/transversion counts:
   `d = -0.5*ln(1-2P-Q) - 0.25*ln(1-2Q)`.
2. FastPhylo's own Python bindings, `fastphylo.distance_matrix(aln,
   model='k2p')` (default `method='ml'`).
3. RapidNJ's own built-in K2P (`rapidnj -i fa -a kim -t d -o m`,
   computed independently, not using FastPhylo at all).
4. `fastdist -D K2P` (default, no extra flags).


(1), (2), and (3) agreed **exactly** on every pair tested (e.g.
`0.079391` for one specific pair). (4) disagreed, and the disagreement
grew with the pair's true distance:

| Pair | Classical/fastphylo-ml/RapidNJ | `fastdist -D K2P` (default) | Ratio |
|---|---|---|---|
| close pair | 0.07939 | 0.08125 | 1.023 |
| ...        | 0.12761 | 0.13221 | 1.036 |
| ...        | 0.19711 | 0.20886 | 1.060 |
| ...        | 0.19539 | 0.20640 | 1.056 |
| more diverged pair | 0.26989 | 0.29258 | **1.084** |

Adding `-N` to the `fastdist -D K2P` call made its output match (1),
(2), and (3) **exactly**, for every pair checked.

## Root cause

`fastdist`'s default path (`no_tstvratio=false`) routes through
`compute_K2P_fixratio()` in `src/fastphylo/dna/Kimura2parameter.cpp`.
That function does a numerical (secant-method) maximum-likelihood
search for the transversion probability `bt_prob`, **given a fixed
transition/transversion ratio `K`** passed in as a parameter - it does
not estimate the ratio from the pair's own data. `K` comes from
`main.cpp`'s `tstvratio` option, which **defaults to `2.0F`**
(`src/c++/apps/fastdist/main.cpp:118`), and is only bypassed (falling
through to the classical closed-form formula at lines 112-129 of
`Kimura2parameter.cpp`) in the degenerate case where a pair's own
transition/transversion ratio is extremely skewed (`a/b > 10000` or
`< 0.00001`).

So for the overwhelming majority of real pairs, `fastdist -D K2P`'s
default behavior is: *assume every pair in the dataset has a
transition/transversion ratio of exactly 2.0, and find the distance
that best explains the observed counts under that assumption* - not
*estimate the distance and the ratio jointly from this pair's own
data*, which is what "K2P distance" means in the standard/textbook
sense, and what `-N`, FastPhylo's own Python `distance_matrix()`, and
RapidNJ's independent implementation all do instead.


When the true (or per-pair-empirical) ts/tv ratio in real data differs
from 2.0 - which is common, and is exactly the situation for
HKY-simulated data, where the true ts/tv-related parameter (kappa) is
whatever the simulator picked, not necessarily 2.0 - this fixed-ratio
assumption biases the ML fit, and the bias compounds with distance
(more divergence -> more informative data about the *true* ratio ->
bigger mismatch against the wrong assumed one -> bigger distance-
estimate correction needed to still explain the observed counts as
well as possible under the wrong ratio).

This is a **design/default-behavior issue**, not an arithmetic bug -
`compute_K2P_fixratio()`'s math is internally consistent for what it's
trying to do (fixed-ratio ML estimation is a legitimate, sometimes-
useful technique, e.g. for stabilizing very short/uninformative
pairs). The problem is that it's silently the *default* for a distance
function whose name and CLI help text ("K2P") give no indication that
anything other than the standard formula is being used, and the `-N`
flag's own help text ("If given fixed ts/tv ratios will not be used")
doesn't explain what happens when it *isn't* given.

## Two separate but related findings, worth deciding on independently

1. **The default itself.** Should `fastdist -D K2P` default to the
   classical per-pair-estimated formula (i.e. today's `-N` behavior),
   with the fixed-ratio ML path (today's default) available as an
   explicit opt-in? This is the behavior change most likely to matter
   to downstream users, since it's the one that silently produces
   biased results for anyone who just runs `fastdist -D K2P` without
   knowing about `-N`.
2. **Documentation.** Regardless of what the default ends up being,
   `-N`'s help text and any `-D K2P`-related docs should explain
   *what changes* - concretely, that the default assumes a fixed ts/tv
   ratio (currently hardcoded to 2.0, overridable via `-T`) rather than
   estimating it per pair, and that this can materially bias distances
   on data whose true ratio differs from that default.



## Suggested fix

Recommended: swap the default for `-D K2P` specifically - make the
per-pair-estimated (classical) computation the default, and require an
explicit flag (e.g. repurpose `-T`/`--tstvratio` as "set this and use
the fixed-ratio ML method instead" rather than "override the assumed
ratio for a method that's already on by default") to opt into the
fixed-ratio path. This matches what a user asking for "K2P" would
reasonably expect, and matches the already-existing classical fallback
inside `compute_K2P_fixratio()` itself (lines 112-129) - the function
already knows how to do the textbook thing, it's just not reaching
that branch except in the degenerate a/b-ratio case.

If changing the default is considered too risky for existing callers
(script/pipeline breakage from a behavior change under the same flags),
the minimum acceptable fix is: update `-N`'s help text (and any
user-facing docs mentioning `-D K2P`) to explicitly state what the
default computes and why it can diverge from the textbook formula, so
this isn't discoverable only by an external user tracing a 2-8%
systematic bias back through three independent reimplementations.

## Verification plan

- Add a regression test comparing `fastdist -D K2P` (whatever the
  chosen default ends up being) against the classical closed-form
  Kimura formula on a range of synthetic pairs spanning low to high
  divergence, asserting agreement within numerical tolerance - this
  specific discrepancy had no test coverage catching it.
- If the default changes: verify `-D JC` and `-D HAMMING` are
  unaffected (both were already confirmed byte-identical with/without
  `-N` on real test data, since neither model has a ts/tv-ratio
  parameter to begin with - this should remain true, just re-confirm
  after any refactor).
- Check whether `TN93` (`-D TN93`) has an analogous issue - it shares
  the `-T`/`--tstvratio`/`-P`/`--pyrtvratio` CLI options (described in
  `--help` as "for the TN model"), suggesting a related or shared code
  path, but this wasn't directly tested as part of this investigation
  (no live caller in the downstream project currently uses `-D TN93` -
  it was tried and rejected there for an unrelated crash on
  near-saturated pairs, so K2P was the only one actually exercised).
  Worth a quick check before considering this done, in case TN93 has
  the same silent-default-bias problem.


## Downstream impact already worked around

The downstream project (`dnctree_testing`) has already applied `-N` to
its own `fastdist -D K2P` calls as a workaround (two call sites, both
named `run_fastdist`, in `experiments/alisim_dna_matrix/Snakefile` and
`experiments/alisim_raxmlgrove_matrix/Snakefile`), so this isn't
blocking their work - this plan is about fixing the root cause in
FastPhylo itself so future callers (including anyone else using
`fastdist -D K2P` without knowing about this) don't need to discover
and work around the same bias independently.