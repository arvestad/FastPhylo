// Measures how sparse the replacement-count matrix N actually is on
// real data, to check whether a sparsity-exploiting rewrite of
// likelihood_slope_curv() (skip full P'(t)/P''(t) rows/cells for
// amino acids that never appear as a "from" state) is worth building.
// Standalone, exploratory. Build:
//   c++ -O2 -std=c++17 -I../include -I/opt/homebrew/Cellar/simde/0.8.2/include \
//     measure_N_sparsity.cpp ../src/fastphylo/protein/ProtSeqCode.cpp \
//     ../src/fastphylo/protein/ProtSeqCompare.cpp -o measure_N_sparsity
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>
#include "fastphylo/protein/ProtSeqCode.hpp"
#include "fastphylo/protein/ProtSeqCompare.hpp"

namespace
{
std::vector<std::string> read_fasta(const std::string &path)
{
    std::ifstream in(path);
    std::vector<std::string> seqs;
    std::string line, current;
    while (std::getline(in, line))
    {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.empty())
            continue;
        if (line[0] == '>')
        {
            if (!current.empty())
            {
                seqs.push_back(current);
                current.clear();
            }
        }
        else
            current += line;
    }
    if (!current.empty())
        seqs.push_back(current);
    return seqs;
}
} // namespace

int main(int argc, char **argv)
{
    std::string path = argc >= 2 ? argv[1] : "../examples/globin_family.fasta";
    std::vector<std::string> seqs = read_fasta(path);
    std::vector<std::vector<std::uint8_t>> encoded(seqs.size());
    for (std::size_t i = 0; i < seqs.size(); i++)
        ProtSeqCode::encode_sequence(seqs[i], encoded[i]);

    std::vector<int> nnz_samples;
    std::vector<int> rows_samples;
    for (std::size_t i = 0; i < seqs.size(); i++)
    {
        for (std::size_t j = i + 1; j < seqs.size(); j++)
        {
            std::vector<std::size_t> tally =
                ProtSeqCode::count_replacement_tally(encoded[i].data(), encoded[i].size(), encoded[j].data(),
                                                       encoded[j].size());
            int nnz = 0;
            std::vector<bool> row_active(20, false);
            for (std::size_t a = 0; a < 20; a++)
            {
                for (std::size_t b = 0; b < 20; b++)
                {
                    if (tally[(a * 20) + b] != 0)
                    {
                        nnz++;
                        row_active[a] = true;
                    }
                }
            }
            int rows = std::count(row_active.begin(), row_active.end(), true);
            nnz_samples.push_back(nnz);
            rows_samples.push_back(rows);
        }
    }

    auto stats = [](std::vector<int> &v, const char *name) {
        std::sort(v.begin(), v.end());
        double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        std::cout << name << ": n=" << v.size() << " min=" << v.front() << " median=" << v[v.size() / 2]
                  << " mean=" << mean << " max=" << v.back() << " p90=" << v[static_cast<std::size_t>(v.size() * 0.9)]
                  << "\n";
    };
    stats(nnz_samples, "nnz(N) out of 400");
    stats(rows_samples, "active rows(N) out of 20");

    return 0;
}
