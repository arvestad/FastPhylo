// Micro-benchmark for Q2 (fastprot ML speedup investigation plan,
// planning/fastprot_ml_speedup_investigation_plan.md): how much does
// MatrixExpm(Q)'s one-time LAPACK decomposition (dgeev_/dgetrf_/
// dgetri_) actually cost, in isolation, across the 5 named models?
// Combine with wall-clock scaling data (bench_decomposition_cost.sh)
// to get decomposition's share of total run time at each dataset
// size - this file only measures the decomposition itself.
//
// Standalone, exploratory - not wired into the CMake build. Build
// manually:
//   c++ -O2 -std=c++17 -I../include \
//     bench_decomposition_cost.cpp \
//     ../src/fastphylo/protein/Matrix.cpp \
//     ../src/fastphylo/protein/ModelMatrix.cpp \
//     -framework Accelerate -o bench_decomposition_cost
//   ./bench_decomposition_cost

#include <algorithm>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include "fastphylo/protein/Matrix.hpp"
#include "fastphylo/protein/ModelMatrix.hpp"

using namespace std::chrono;

namespace
{

struct Stats
{
    double median_us;
    double q1_us;
    double q3_us;
};

Stats summarize(std::vector<double> &samples_us)
{
    std::sort(samples_us.begin(), samples_us.end());
    std::size_t n = samples_us.size();
    Stats s;
    s.median_us = samples_us[n / 2];
    s.q1_us = samples_us[n / 4];
    s.q3_us = samples_us[(3 * n) / 4];
    return s;
}

} // namespace

int main()
{
    const int WARMUP = 3;
    const int REPS = 51; // decomposition is a bit noisier than the pure-multiply Q1 benchmark - more reps

    struct
    {
        model_type type;
        const char *name;
    } models[] = {{wag, "wag"}, {jtt, "jtt"}, {day, "day"}, {mvr, "mvr"}, {lg, "lg"}};

    std::cout << "model,median_us,q1_us,q3_us\n";
    for (auto &m : models)
    {
        Matrix Q = get_model_matrix(m.type);

        for (int i = 0; i < WARMUP; i++)
        {
            MatrixExpm decomp(Q);
            (void)decomp;
        }

        std::vector<double> samples;
        samples.reserve(REPS);
        for (int i = 0; i < REPS; i++)
        {
            auto t0 = high_resolution_clock::now();
            MatrixExpm decomp(Q);
            auto t1 = high_resolution_clock::now();
            samples.push_back(duration<double, std::micro>(t1 - t0).count());
            (void)decomp;
        }
        Stats s = summarize(samples);
        std::cout << m.name << "," << s.median_us << "," << s.q1_us << "," << s.q3_us << "\n";
    }

    return 0;
}
