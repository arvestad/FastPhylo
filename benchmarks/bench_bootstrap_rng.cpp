// Measures candidate replacements for bootstrap_sequences()'s
// `seqlen * 1.0 * rand() / (RAND_MAX + 1.0)` position sampling, after
// the first attempt (mt19937 + std::uniform_int_distribution)
// measured SLOWER end-to-end (0.41s -> 0.57s on a -b 50000 stress
// test) than the rand() baseline it was meant to replace -
// uniform_int_distribution's bias-avoidance machinery costs more than
// rand()'s external-call overhead saved. Isolates each candidate's
// per-call cost directly rather than guessing from the earlier,
// contradicted assumption that any <random> engine trivially wins.
//
// Result adopted: mt19937 + the same float-scaling trick the old
// rand()-based code already used (skip uniform_int_distribution
// entirely) - 3.71ns/call vs rand()'s 7.18ns/call, ~1.94x faster,
// standard-library only (fastphylo/core/random_utils.cpp). xorshift32
// measured faster still (2.04ns/call) but was deliberately not chosen -
// its base 3-shift form is a known-linear generator (fails binary-rank/
// linear-complexity tests, per Marsaglia's own later xorshift+ work and
// Vigna's xoshiro papers, both of which add a scrambling step
// specifically to fix this) - a standard, well-vetted engine was
// preferred over a faster but non-standard, weaker one.
//
// Standalone, exploratory - not wired into the CMake build. Build:
//   c++ -O2 -std=c++17 bench_bootstrap_rng.cpp -o bench_bootstrap_rng
//   ./bench_bootstrap_rng

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

using namespace std::chrono;

namespace
{
constexpr int N = 300; // matches ml_bench_2x300.fasta's seqlen
constexpr long REPS = 50'000'000;

double median_of(std::vector<double> v)
{
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

template <typename F> double time_it(const char *name, F &&f)
{
    volatile long sink = 0;
    // warmup
    for (long i = 0; i < 1'000'000; i++)
    {
        sink += f();
    }
    std::vector<double> samples;
    for (int rep = 0; rep < 5; rep++)
    {
        auto t0 = high_resolution_clock::now();
        for (long i = 0; i < REPS; i++)
        {
            sink += f();
        }
        auto t1 = high_resolution_clock::now();
        samples.push_back(duration<double, std::milli>(t1 - t0).count());
    }
    double ms = median_of(samples);
    std::cout << name << ": " << ms << " ms for " << REPS << " calls (" << (ms * 1e6 / REPS) << " ns/call)"
              << " [sink=" << sink << "]\n";
    return ms;
}

// xorshift32 - minimal state (4 bytes), no tempering step, fully
// inlinable. Not cryptographic/high-quality, but bootstrap resampling
// doesn't need that - it needs a statistically reasonable, fast
// source of uniform integers.
struct Xorshift32
{
    std::uint32_t state;
    std::uint32_t operator()()
    {
        std::uint32_t x = state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        state = x;
        return x;
    }
};

} // namespace

int main()
{
    srand(1);
    time_it("rand() + float scaling (current baseline)",
            [&]() { return static_cast<int>(N * 1.0 * rand() / (RAND_MAX + 1.0)); });

    std::mt19937 mt(1);
    time_it("mt19937 + uniform_int_distribution (first attempt)", [&]() {
        std::uniform_int_distribution<int> dist(0, N - 1);
        return dist(mt);
    });

    std::mt19937 mt2(1);
    time_it("mt19937 + float scaling (skip uniform_int_distribution)",
            [&]() { return static_cast<int>(N * 1.0 * mt2() / (static_cast<double>(mt2.max()) + 1.0)); });

    Xorshift32 xs{1};
    time_it("xorshift32 + float scaling",
            [&]() { return static_cast<int>(N * 1.0 * xs() / (static_cast<double>(UINT32_MAX) + 1.0)); });

    std::minstd_rand lcg(1);
    time_it("std::minstd_rand + float scaling",
            [&]() { return static_cast<int>(N * 1.0 * lcg() / (static_cast<double>(lcg.max()) + 1.0)); });

    return 0;
}
