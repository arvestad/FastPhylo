// Dumps every named model's Q matrix and equilibrium frequencies as
// plain numbers, for cross-checking against IQ-TREE2's embedded model
// data (model/modelprotein.cpp) via compare_to_iqtree2.py - written
// after finding JTT's Q(0,0) was hardcoded 100x too large (see
// tests/ModelMatrix_test.cpp and planning/fastprot_ml_speedup_
// implementation_plan.md). One line per model: "MODELNAME_Q" then 400
// row-major values, then "MODELNAME_EQ" then 20 values.
//
// Standalone, exploratory - not wired into the CMake build. Build:
//   c++ -O2 -std=c++17 -I../include -I<path-to-eigen3-include-dir> \
//     dump_model_matrices.cpp ../src/fastphylo/protein/Matrix.cpp \
//     ../src/fastphylo/protein/ModelMatrix.cpp -framework Accelerate \
//     -o dump_model_matrices   # drop -framework Accelerate, add
//                               # -lblas -llapack on Linux
//   ./dump_model_matrices > /tmp/fastprot_matrices.txt
//
// Then, to re-run the comparison:
//   curl -sL https://raw.githubusercontent.com/iqtree/iqtree2/master/model/modelprotein.cpp \
//     -o /tmp/iqtree_modelprotein.cpp
//   python3 compare_to_iqtree2.py
#include <iomanip>
#include <iostream>
#include "fastphylo/protein/Matrix.hpp"
#include "fastphylo/protein/ModelMatrix.hpp"

int main()
{
    struct
    {
        model_type type;
        const char *name;
    } models[] = {{wag, "WAG"}, {jtt, "JTT"}, {day, "DAYHOFF"}, {mvr, "MVR"}, {lg, "LG"}};

    std::cout << std::setprecision(10);
    for (auto &m : models)
    {
        Matrix Q = get_model_matrix(m.type);
        DblVec eq = get_model_vec(m.type);
        std::cout << m.name << "_Q\n";
        for (std::size_t row = 0; row < 20; row++)
        {
            for (std::size_t col = 0; col < 20; col++)
            {
                std::cout << Q(row, col) << (col + 1 < 20 ? " " : "\n");
            }
        }
        std::cout << m.name << "_EQ\n";
        for (std::size_t i = 0; i < eq.size(); i++)
        {
            std::cout << eq[i] << (i + 1 < eq.size() ? " " : "\n");
        }
    }
    return 0;
}
