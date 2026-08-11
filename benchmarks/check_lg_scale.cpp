// Checks whether LG's embedded rate matrix (ModelMatrix.cpp) is on
// the same raw magnitude scale as WAG/JTT/Dayhoff/MVR's, or a
// different one - prompted by analyze_convergence.cpp finding LG
// converges very differently from the other 4 models after the
// PAM_TO_SUBSTITUTIONS_PER_SITE rescale. This off-diagonal-sum check
// gave the first clear signal (LG=19.69 vs ~0.19 for the others) but
// the actual fix (mean_substitution_rate(), ModelMatrix.hpp/.cpp,
// 2026-08-11) uses the precise, standard "mean substitution rate"
// definition instead - kept here as the exploratory evidence that
// found the bug, not as the reference computation.
#include <iostream>
#include <numeric>
#include <vector>
#include "fastphylo/protein/Matrix.hpp"
#include "fastphylo/protein/ModelMatrix.hpp"

int main()
{
    struct
    {
        model_type type;
        const char *name;
    } models[] = {{wag, "WAG"}, {jtt, "JTT"}, {day, "Dayhoff"}, {mvr, "MVR"}, {lg, "LG"}};

    for (auto &m : models)
    {
        Matrix Q = get_model_matrix(m.type);
        double off_diag_sum = 0;
        double max_off_diag = 0;
        for (std::size_t r = 0; r < 20; r++)
        {
            for (std::size_t c = 0; c < 20; c++)
            {
                if (r != c)
                {
                    off_diag_sum += std::fabs(Q(r, c));
                    max_off_diag = std::max(max_off_diag, std::fabs(Q(r, c)));
                }
            }
        }
        std::cout << m.name << ": off_diag_sum=" << off_diag_sum << " max_off_diag_entry=" << max_off_diag << "\n";
    }
    return 0;
}
