#pragma once

#include <cstddef>

//! Seeds the process-wide generator uniform_random_index() draws from -
//! call once, before generating any bootstrap replicate. Replaces the
//! old srand() call: rand()/srand() are deprecated by the C++ Core
//! Guidelines (SL.random.1) and, in practice here, were measured as a
//! real cost - a non-inlinable external libc call with implicit global
//! state, rather than the header-only, inlinable engines <random>
//! offers (planning/fastprot_ml_speedup_implementation_plan.md's
//! profiling round found it ~21% of wall time in a bootstrap-heavy
//! workload). Kept as simple global state (mirroring rand()/srand()'s
//! existing shape) rather than an explicitly-threaded generator object,
//! since bootstrap resampling is inherently sequential/single-threaded
//! today and this is a direct swap-in, not a new architecture.
void seed_rng(unsigned int seed);

//! A uniformly-distributed integer in [0, n) - replaces the
//! `n * 1.0 * rand() / (RAND_MAX + 1.0)` scaling pattern duplicated
//! across ProtSeqUtils.cpp/Sequence.cpp/Sequences2DistanceMatrix.cpp's
//! bootstrap-resampling functions (all sample a position in
//! [0, seqlen) this same way). std::uniform_int_distribution also
//! avoids any modulo-bias concern by construction, not just by the
//! float-division trick the old pattern relied on.
int uniform_random_index(std::size_t n);
