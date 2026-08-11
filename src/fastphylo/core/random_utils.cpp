#include "fastphylo/core/random_utils.hpp"
#include <random>

namespace
{
std::mt19937 g_rng;
}

void seed_rng(unsigned int seed)
{
    g_rng.seed(seed);
}

int uniform_random_index(std::size_t n)
{
    // Not std::uniform_int_distribution: measured 2.45x SLOWER than
    // the old rand()-based code it was meant to replace (17.6ns/call
    // vs 7.2ns/call) - its bias-avoidance machinery costs more than
    // rand()'s external-call overhead saved. This is the same float-
    // scaling approach the old rand()-based code used (n * 1.0 *
    // draw() / (max()+1.0)), just fed by mt19937 instead of rand() -
    // measured 3.71ns/call, 1.94x faster than the original rand()
    // code (benchmarks/bench_bootstrap_rng.cpp). Reintroduces the same
    // theoretical, negligible-at-this-n modulo-style bias the old code
    // already had (n=~300 vs mt19937's 2^32 range) rather than paying
    // for uniform_int_distribution's stricter unbiasing.
    return static_cast<int>(static_cast<double>(n) * g_rng() / (static_cast<double>(g_rng.max()) + 1.0));
}
