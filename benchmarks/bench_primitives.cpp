// Micro-benchmark: old (toupper()+compare / getAAInd()-based) vs new
// (ProtSeqCode-encoded SIMD) comparison primitives, in isolation.
// Phase 5 of the speed2026a protein-distance speedup (plan.md).
//
// Methodology: for each sequence length, generate one fixed-seed pair at
// a representative divergence level, run 3 warm-up calls then 21 timed
// calls, report the median and interquartile range (not a single
// timing, to avoid being misled by system noise - plan.md's Phase 5
// requirement). Encoding (ProtSeqCode::encode_sequence) happens once per
// sequence, outside the timed region, matching how Phase 3 wired the new
// primitive into fastprot (encode once per sequence, not per pair).
//
// Build via CMake (target bench_primitives) or manually:
//   c++ -O2 -std=c++17 -I../include -I/opt/homebrew/include \
//     bench_primitives.cpp ../src/fastphylo/protein/ProtSeqCode.cpp \
//     ../src/fastphylo/protein/ProtSeqCompare.cpp -o bench_primitives
//   ./bench_primitives

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include "fastphylo/protein/ProtSeqCode.hpp"
#include "fastphylo/protein/ProtSeqCompare.hpp"

using namespace ProtSeqCode;
using namespace std::chrono;

namespace {

constexpr std::string_view CANONICAL = "ARNDCQEGHILKMFPSTWYV";

std::string random_seq(std::size_t len, unsigned int *seed) {
  std::string s;
  s.reserve(len);
  for (std::size_t i = 0; i < len; i++)
    s += CANONICAL[rand_r(seed) % CANONICAL.size()];
  return s;
}

std::string diverge(const std::string &s, double rate, unsigned int *seed) {
  std::string out = s;
  for (std::size_t i = 0; i < out.size(); i++)
    if ((rand_r(seed) / (double)RAND_MAX) < rate)
      out[i] = CANONICAL[rand_r(seed) % CANONICAL.size()];
  return out;
}

// Mirrors count_id_dist() (ProtSeqUtils.cpp) exactly.
double old_count_id_dist(const std::string &s1, const std::string &s2) {
  double id = 0;
  for (std::size_t i = 0; i < s1.size(); i++)
    if (std::toupper(static_cast<unsigned char>(s1[i])) == std::toupper(static_cast<unsigned char>(s2[i])))
      id++;
  return id / static_cast<double>(s1.size());
}

std::size_t old_getAAInd(char c) {
  std::string::size_type pos = CANONICAL.find(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
  return pos == std::string::npos ? 100 : pos;
}

// Mirrors count_replacements() (ProtSeqUtils.cpp) exactly.
std::vector<std::size_t> old_count_replacements(const std::string &s1, const std::string &s2) {
  std::vector<std::size_t> counts(std::size_t{20} * 20, 0);
  for (std::size_t i = 0; i < s1.size(); i++) {
    std::size_t c1 = old_getAAInd(s1[i]);
    std::size_t c2 = old_getAAInd(s2[i]);
    if (c1 != 100 && c2 != 100)
      counts[c1 * 20 + c2]++;
  }
  return counts;
}

struct Stats {
  double median_us;
  double q1_us;
  double q3_us;
};

Stats summarize(std::vector<double> &samples_us) {
  std::sort(samples_us.begin(), samples_us.end());
  std::size_t n = samples_us.size();
  Stats s;
  s.median_us = samples_us[n / 2];
  s.q1_us = samples_us[n / 4];
  s.q3_us = samples_us[(3 * n) / 4];
  return s;
}

template <typename Fn>
Stats time_calls(Fn fn, int warmup, int reps) {
  for (int i = 0; i < warmup; i++)
    fn();
  std::vector<double> samples;
  samples.reserve(reps);
  for (int i = 0; i < reps; i++) {
    auto t0 = high_resolution_clock::now();
    fn();
    auto t1 = high_resolution_clock::now();
    samples.push_back(duration<double, std::micro>(t1 - t0).count());
  }
  return summarize(samples);
}

} // anonymous namespace

int main() { // NOLINT(bugprone-exception-escape) - a benchmark binary crashing on an unexpected exception is a fine, loud failure mode.
  const int WARMUP = 3;
  const int REPS = 21;
  const double DIVERGENCE = 0.25; // representative mid-range divergence
  const std::size_t lengths[] = { 50, 300, 2000 }; // short / typical / long

  std::cout << "primitive,length,old_median_us,old_q1_us,old_q3_us,new_median_us,new_q1_us,new_q3_us,speedup\n";

  for (std::size_t li = 0; li < sizeof(lengths) / sizeof(lengths[0]); li++) {
    std::size_t len = lengths[li];
    unsigned int seed = 42 + static_cast<unsigned int>(len);
    std::string s1 = random_seq(len, &seed);
    std::string s2 = diverge(s1, DIVERGENCE, &seed);

    std::vector<std::uint8_t> c1, c2;
    encode_sequence(s1, c1);
    encode_sequence(s2, c2);

    volatile double sink_d = 0;
    Stats old_id = time_calls([&]() { sink_d += old_count_id_dist(s1, s2); }, WARMUP, REPS);
    Stats new_id = time_calls([&]() {
      sink_d += count_id_fraction(c1.data(), c1.size(), c2.data(), c2.size());
    }, WARMUP, REPS);
    std::cout << "count_id_dist," << len << ","
              << old_id.median_us << "," << old_id.q1_us << "," << old_id.q3_us << ","
              << new_id.median_us << "," << new_id.q1_us << "," << new_id.q3_us << ","
              << (old_id.median_us / new_id.median_us) << "\n";

    volatile std::size_t sink_s = 0;
    Stats old_rep = time_calls([&]() {
      std::vector<std::size_t> t = old_count_replacements(s1, s2);
      sink_s += t[0];
    }, WARMUP, REPS);
    Stats new_rep = time_calls([&]() {
      std::vector<std::size_t> t = count_replacement_tally(c1.data(), c1.size(), c2.data(), c2.size());
      sink_s += t[0];
    }, WARMUP, REPS);
    std::cout << "count_replacements," << len << ","
              << old_rep.median_us << "," << old_rep.q1_us << "," << old_rep.q3_us << ","
              << new_rep.median_us << "," << new_rep.q1_us << "," << new_rep.q3_us << ","
              << (old_rep.median_us / new_rep.median_us) << "\n";

    if (sink_d == 123456.0 && sink_s == 123456) std::cerr << "unreachable\n"; // prevent over-optimization
  }

  return 0;
}
