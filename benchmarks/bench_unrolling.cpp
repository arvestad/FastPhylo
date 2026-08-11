// Tests whether raising EIGEN_UNROLLING_LIMIT (default 110, checked
// against the installed Eigen 3.4.0 headers - a 20x20x20 dense
// multiply is ~8000 multiply-adds, well past that default) changes
// measured performance for the exact operations likelihood_slope_curv()
// (MaximumLikelihood.cpp) performs every Newton iteration.
//
// EIGEN_UNROLLING_LIMIT must be defined before Eigen is first
// #included in a translation unit, and mixing different values of it
// across TUs linked into one binary would be an ODR violation (the
// same template instantiation would have two different bodies) - so
// this is built as two entirely separate standalone binaries (one per
// -DEIGEN_UNROLLING_LIMIT value), not two code paths in one binary.
//
// Result (2026-08-11): no benefit. The isolated 20x20 multiply's
// median time was identical (0.333us) at both the default limit (110)
// and a much higher one (20000, comfortably above the ~8000
// multiply-adds a 20x20x20 product costs) - raising the macro didn't
// change the generated code for this operation at all. The full-
// iteration numbers were, if anything, marginally worse under the
// raised limit, most plausibly just inter-process measurement noise
// (CPU frequency/cache state differing between two separate binary
// invocations) rather than a real regression. Not adopted.
//
// Standalone, exploratory - not wired into the CMake build. Build
// (twice, with different -DEIGEN_UNROLLING_LIMIT):
//   c++ -O2 -std=c++17 -I<path-to-eigen3-include-dir> \
//     -DEIGEN_UNROLLING_LIMIT=20000 \
//     bench_unrolling.cpp -o bench_unrolling_unrolled
//   c++ -O2 -std=c++17 -I<path-to-eigen3-include-dir> \
//     bench_unrolling.cpp -o bench_unrolling_default
//   ./bench_unrolling_default && ./bench_unrolling_unrolled

#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>
#include <Eigen/Dense>

using namespace std::chrono;

namespace
{

struct Stats
{
    double median_us;
};

Stats summarize(std::vector<double> &samples_us)
{
    std::sort(samples_us.begin(), samples_us.end());
    return {samples_us[samples_us.size() / 2]};
}

template <typename Fn> Stats time_calls(Fn fn, int warmup, int reps)
{
    for (int i = 0; i < warmup; i++)
    {
        fn();
    }
    std::vector<double> samples;
    samples.reserve(reps);
    for (int i = 0; i < reps; i++)
    {
        auto t0 = high_resolution_clock::now();
        fn();
        auto t1 = high_resolution_clock::now();
        samples.push_back(duration<double, std::micro>(t1 - t0).count());
    }
    return summarize(samples);
}

} // namespace

int main()
{
    const int WARMUP = 5;
    const int REPS = 51;

    // Representative fixed matrices - values don't matter for timing,
    // just need to be non-degenerate.
    Eigen::Matrix<float, 20, 20> eigenvectors = Eigen::Matrix<float, 20, 20>::Random();
    Eigen::Matrix<float, 20, 20> eigenvectors_inv = Eigen::Matrix<float, 20, 20>::Random();
    Eigen::Matrix<float, 20, 1> eigenvalues = -Eigen::Matrix<float, 20, 1>::Random().cwiseAbs();
    Eigen::Matrix<float, 20, 20> Q = Eigen::Matrix<float, 20, 20>::Random();
    Eigen::Matrix<float, 20, 20> Q2 = Q * Q;
    Eigen::Matrix<float, 20, 20> N = Eigen::Matrix<float, 20, 20>::Random().cwiseAbs();
    float t = 0.35f;

    volatile double sink = 0;

    Stats full_iteration = time_calls(
        [&]() {
            // Mirrors likelihood_slope_curv() exactly.
            Eigen::Matrix<float, 20, 1> scale = (eigenvalues * t).array().exp();
            Eigen::Matrix<float, 20, 20> pt = eigenvectors * scale.asDiagonal() * eigenvectors_inv;
            Eigen::Matrix<float, 20, 20> p1t = pt * Q;
            Eigen::Matrix<float, 20, 20> p2t = pt * Q2;

            Eigen::Array<float, 20, 20> p = pt.array().max(1e-30f);
            Eigen::Array<float, 20, 20> ratio = p1t.array() / p;
            Eigen::Array<double, 20, 20> slope_terms = (N.array() * ratio).cast<double>();
            Eigen::Array<double, 20, 20> curv_terms = (N.array() * (p2t.array() / p - ratio.square())).cast<double>();
            sink += slope_terms.sum() + curv_terms.sum();
        },
        WARMUP, REPS);

    Stats mult_alone = time_calls(
        [&]() {
            Eigen::Matrix<float, 20, 20> r = eigenvectors * eigenvectors_inv;
            sink += r(0, 0);
        },
        WARMUP, REPS);

#ifdef EIGEN_UNROLLING_LIMIT
    std::cout << "EIGEN_UNROLLING_LIMIT=" << EIGEN_UNROLLING_LIMIT << "\n";
#else
    std::cout << "EIGEN_UNROLLING_LIMIT=<default>\n";
#endif
    std::cout << "full_likelihood_slope_curv_iteration_median_us," << full_iteration.median_us << "\n";
    std::cout << "single_20x20_multiply_median_us," << mult_alone.median_us << "\n";

    if (sink == 123456.0)
    {
        std::cerr << "unreachable\n";
    }
    return 0;
}
