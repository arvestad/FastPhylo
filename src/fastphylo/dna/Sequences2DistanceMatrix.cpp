//--------------------------------------------------
//
// File: Sequences2DistanceMatrix.cpp
//
// Author: Isaac Elias,(addition and changes by Mehmood Alam Khan)
// e-mail: isaac@nada.kth.se, malagori@kth.se
//
// cvs: $Id: Sequences2DistanceMatrix.cpp,v 1.51 2006/12/31 11:17:52 isaac Exp $
//
//--------------------------------------------------

#include "fastphylo/dna/Sequences2DistanceMatrix.hpp"
#include "fastphylo/core/Exception.hpp"
#include "fastphylo/dna/DNA_b128_String.hpp"
#include "fastphylo/core/phylip_interleaved_reader.hpp"
#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "fastphylo/core/stl_utils.hpp"
#include "fastphylo/core/log_utils.hpp"
#include "fastphylo/core/file_utils.hpp"
#include <cfloat>
#include <cstdio>

using namespace std;

void Sequences2DNA_b128(std::vector<Sequence> &seqs, std::vector<DNA_b128_String> &b128)
{

    b128.resize(seqs.size());
    for (size_t i = 0; i < seqs.size(); i++)
    {
        size_t cap = seqs[i].seq.length() + 1;
        b128[i].reInitiate(cap);
        b128[i].append(seqs[i].seq);
    }
}

void DNA_b128_StringsFromPHYLIP(istream &fin, std::vector<std::string> &names,
                                std::vector<DNA_b128_String> &b128_strings)
{

    int numSequences;
    int seqlen;
    const int MAXLINE = 16384;
    std::array<char, MAXLINE> line{};
    do
    { // skip lines that does not contain two integers
        fin.getline(line.data(), MAXLINE);
        // NOLINTNEXTLINE(bugprone-unchecked-string-to-number-conversion) - the while condition below is exactly that
        // check.
    } while (sscanf(line.data(), "%d %d", &numSequences, &seqlen) != 2);

    names.clear();
    names.reserve(numSequences);
    b128_strings.resize(numSequences);
    for (int i = 0; i < numSequences; i++)
    {
        b128_strings[i].reInitiate(seqlen);
        names.emplace_back();
    }
    // phylip has name lenght 10.
    // read the names and map the sequences onto the tree
    std::array<char, 11> tmpName{};
    for (int i = 0; i < numSequences; i++)
    {
        DNA_b128_String &s = b128_strings[i];

        fin.getline(tmpName.data(), 11); // reads atmost 10 chars
        // skip lines without 10 chars per line
        if (!fin.fail())
        {
            if (fin.eof())
                THROW_EXCEPTION("Unexpected reading format fin.eof() == " << fin.eof());
            i--;
            continue;
        }
        appendUntil(names[i], tmpName.data(), fin.gcount(), ' ');

        fin.clear();
        fin.getline(line.data(), MAXLINE);

        while (fin.fail() && fin.gcount() == MAXLINE - 1)
        { // didn't read all the line
            s.append(line.data());
            fin.clear();
            fin.getline(line.data(), MAXLINE);
        }
        if (!fin.fail())
        { // we read it all including the newline char unless it ended with eof
            s.append(line.data());
        }
        else // fail
            THROW_EXCEPTION("Unexpected reading format fin.eof() == " << fin.eof());
    } // end for loop

    // Mehmood changes starts here
    // The sequences aren't neccesarily on one line but my be spread out interleaving
    // over several lines. Therefore we read until seqlen chars have been read.
    phylipReadInterleavedContinuation(
        fin, numSequences, static_cast<size_t>(seqlen),
        [&b128_strings](int i) { return static_cast<size_t>(b128_strings[i].getNumChars()); },
        [&](istream &in, int i) {
            DNA_b128_String &s = b128_strings[i];
            while (true)
            {
                in.getline(line.data(), MAXLINE);
                std::string myStr = line.data();
                if (myStr.empty())
                {
                    if (in.eof())
                    {
                        THROW_EXCEPTION("Sequence not of correct length: " << names[i] << "    length is "
                                                                           << b128_strings[i].getNumChars());
                    }
                    continue;
                }
                if (!in.fail() || in.gcount() != MAXLINE - 1)
                {
                    s.append(line.data());
                }
                break;
            }
        });

    // CHECK THAT ALL STRINGS HAVE THE SAME LENGTH
    for (int i = 0; i < numSequences; i++)
    {
        if (b128_strings[i].getNumChars() != seqlen)
        {
            THROW_EXCEPTION("Sequence not of correct length: " << names[i] << "    length is "
                                                               << b128_strings[i].getNumChars());
            DNA_b128_String &s = b128_strings[i];
            cerr << s << endl;
        }
    }
    //-----------------------------
}
//---------------------------------------------------------
void bootstrapSequences(const std::vector<Sequence> &seqs, std::vector<DNA_b128_String> &b128_strings)
{

    // ensure capacity in bootsequences
    b128_strings.resize(seqs.size());

    const size_t seqlen = seqs[0].seq.length();
    for (size_t i = 0; i < seqs.size(); i++)
    {
        b128_strings[i].reInitiate(seqlen);
    }

    // Do the bootstrapping
    size_t pos = 0;
    vector<int> samplePositions(seqlen);
    const size_t BUFFSIZE = (16383 > seqlen ? seqlen : 16383); // 2^14=16384
    std::vector<char> buff(BUFFSIZE + 1);
    buff[BUFFSIZE] = '\0';

    // Was two near-identical branches keyed on `32 < seqlen`: one doing
    // this same fill in stride-32 chunks (an apparent unrolling
    // attempt that, since rand() is an inherently serial call, visits
    // positions 0..seqlen-1 in the same order and calls rand() the
    // same number of times as the plain loop below - chunking it
    // changed nothing observable), the other this plain loop, used
    // only for seqlen<=32. Proven equivalent and unified into one path.
    for (pos = 0; pos < seqlen; pos++)
    {
        samplePositions[pos] = static_cast<int>(seqlen * 1.0 * rand() / (RAND_MAX + 1.0));
    }

    // Same merge for the sampled-sequence-copy step below: BUFFSIZE is
    // already capped at seqlen (min(16383,seqlen)), so for seqlen<=32
    // (where the "stride<seqlen" branch used to skip this chunking
    // entirely) the chunk loop below runs zero iterations and the tail
    // loop alone copies the whole sequence in one append() - exactly
    // what the removed branch did by hand.
    for (size_t seq = 0; seq < seqs.size(); seq++)
    {
        const string &s = seqs[seq].seq;
        pos = 0;
        for (; pos < seqlen - BUFFSIZE; pos += BUFFSIZE)
        {
            for (size_t i = 0; i < BUFFSIZE; i++)
            {
                buff[i] = s[samplePositions[pos + i]];
            }
            b128_strings[seq].append(buff.data());
        }
        size_t i;
        for (i = 0; pos < seqlen; i++, pos++)
        {
            buff[i] = s[samplePositions[pos]];
        }
        buff[i] = '\0';
        b128_strings[seq].append(buff.data());
    }
}

//----------------------------------------------------------

void fillMatrix(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs, sequence_translation_model trans_model)
{

    const size_t numSequences = seqs.size();
    dm.resize(numSequences);

    // base frequences
    DNA_b128_String::base_frequences freqs = seqs[0].getBaseFrequences();
    for (size_t i = 1; i < numSequences; i++)
    {
        DNA_b128_String::base_frequences tmpfreqs = seqs[i].getBaseFrequences();
        freqs.num_As_ += tmpfreqs.num_As_;
        freqs.num_Cs_ += tmpfreqs.num_Cs_;
        freqs.num_Gs_ += tmpfreqs.num_Gs_;
        freqs.num_Ts_ += tmpfreqs.num_Ts_;
        freqs.num_unknowns_ += tmpfreqs.num_unknowns_;
        freqs.num_ambiguities_ += tmpfreqs.num_ambiguities_;
    }

    // compute ambiguities probabilities
    if (!trans_model.no_ambiguities)
    {
        if (trans_model.use_base_freqs)
        {
            for (size_t i = 0; i < numSequences; i++)
            {
                seqs[i].calcAmbiguityProbabilities(freqs.num_As_, freqs.num_Cs_, freqs.num_Gs_, freqs.num_Ts_);
            }
        }
        else
        {
            for (size_t i = 0; i < numSequences; i++)
            {
                seqs[i].calcAmbiguityProbabilitiesUNIFORM();
            }
        }
    }

    // CALL DISTANCE COMPUTATION
    switch (trans_model.model)
    {
    case HAMMING_DISTANCE:
        fillMatrix_Hamming(dm, seqs, trans_model);
        break;
    case JC:
        fillMatrix_JC(dm, seqs, trans_model);
        break;
    case K2P:
        fillMatrix_K2P(dm, seqs, trans_model);
        break;
    case TN93:
        fillMatrix_TN93(dm, seqs, freqs, trans_model);
        break;
    default:
        THROW_EXCEPTION("Non handled model: " << trans_model.model);
    }
}

// distance_matrix_refactor_plan.md, Phase C: fillMatrix_{Hamming,JC,
// K2P,TN93} and fillMatrixRow_{Hamming,JC,K2P,TN93} (8 functions
// total, ~90% identical to each other) are consolidated below into
// three pieces of code instead:
//
// - fillMatrixML<RawDistance>(): shared by fillMatrix_JC/_K2P/_TN93.
//   JC/K2P/TN93 share one skeleton exactly (per-pair raw distance ->
//   ML distance, closest-neighbor tracking, an ambiguity-resolution
//   pass using that ML distance, then a correction pass) - only the
//   raw-distance type (simple_string_distance for JC/K2P,
//   TN_string_distance for TN93) and the three model-specific
//   callables (computeRaw/computeML/correctAmbiguities) differ.
// - fillMatrix_Hamming() stays its own function: its ambiguity
//   handling is genuinely different in shape, not just formula -
//   compute_Hamming_distance() returns a plain float (not an
//   ML_string_distance), its closest-neighbor resolution calls the
//   simpler resolveAmbiguities() (not
//   resolveAmbiguitiesUsingTransitionProbabilities()), and its
//   correction pass borrows compute_JC() as an ML estimate regardless
//   of the model. Forcing this into the same template as the other
//   three would need enough extra parameterization (a
//   "closest-neighbor resolver" callable, a separate
//   "correction-estimate" callable, a storage-shape switch) that it
//   would be harder to read than the duplication it replaced -
//   judged not worth it.
// - fillMatrixRowGeneric<RawDistance>(): shared by all four
//   fillMatrixRow_* wrappers, including Hamming. Unlike the
//   full-matrix versions, none of the four *row*-streaming originals
//   had a live ambiguity-resolution pass - all four had it entirely
//   commented out already (dead code, deleted here rather than kept
//   as an inert comment four times over) - so once that's gone, even
//   Hamming's row version reduces to exactly the same shape as the
//   other three (raw distance -> a plain float, no ML type or
//   closest-neighbor bookkeeping needed at all), letting all four
//   unify into one template.
//
// Bug found and fixed during this consolidation: fillMatrixRow_JC's
// `mem_eff_flag` branch was inverted relative to its three siblings
// (`if (mem_eff_flag)` instead of `if (!mem_eff_flag)`). Verified
// live (not just read): `fastdist -D JC -e -O phylip` on an 8-sequence
// test set printed zeros for the entire lower triangle instead of the
// real distances, while the same command with `-D K2P` produced
// output byte-identical to the non-memory-efficient full matrix (the
// correct reference). fillMatrixRowGeneric() below uses the
// three-out-of-four (correct) direction uniformly, so this bug cannot
// recur in any of the four models.
//
// The `row = -1` (relying on unsigned wraparound to make `row+1`
// become `0`) trick every original fillMatrixRow_* used to switch
// between "start at row+1" and "start at 0" is replaced by an
// explicit `startCol` variable below - a reader doesn't have to know
// unsigned integer wraparound is happening on purpose to understand
// the control flow, and it can't silently invert again the way JC's
// copy did.

namespace
{

// Shared by fillMatrix_JC/_K2P/_TN93 - see the consolidation comment
// above for why fillMatrix_Hamming isn't included here.
template <class RawDistance, class ComputeRawFn, class ComputeMLFn, class CorrectAmbiguitiesFn>
// NOLINTNEXTLINE(readability-function-cognitive-complexity) - same shape as the fillMatrix_{JC,K2P,TN93} it replaces,
// all three were already at this same complexity before consolidation.
void fillMatrixML(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs, const sequence_translation_model &trans_model,
                  ComputeRawFn computeRaw, ComputeMLFn computeML, CorrectAmbiguitiesFn correctAmbiguities)
{
    const size_t numSequences = seqs.size();
    dm.resize(numSequences);

    // extended information used to compute ambiguities
    static DistanceMatrix<int, pair<RawDistance, ML_string_distance>, empty_Data_init<int>, empty_Data_printOn<int>,
                          empty_Data_init<pair<RawDistance, ML_string_distance>>,
                          empty_Data_printOn<pair<RawDistance, ML_string_distance>>>
        extendedDistanceInfo(numSequences);
    extendedDistanceInfo.resize(numSequences);

    for (size_t i = 0; i < numSequences; i++)
    {
        dm.setDistance(i, i, 0);
        size_t closestNeig = i;
        float closestDist = FLT_MAX;

        DNA_b128_String &si = seqs[i];

        // if has ambig find closest string and resolve.
        if (si.hasAmbiguities())
        {
            for (size_t k = 0; k < i; k++)
            {
                float dist = dm.getDistance(k, i);
                if (dist < closestDist && dist >= 0)
                {
                    closestDist = dist;
                    closestNeig = k;
                }
            }
        }

        for (size_t j = i + 1; j < numSequences; j++)
        {
            RawDistance sd = computeRaw(si, seqs[j]);
            ML_string_distance ml_dist = computeML(sd);

            dm.setDistance(i, j, ml_dist.distance);
            extendedDistanceInfo.setDistance(i, j, pair<RawDistance, ML_string_distance>(sd, ml_dist));

            if (ml_dist.distance < closestDist && ml_dist.distance >= 0)
            {
                closestDist = ml_dist.distance;
                closestNeig = j;
            }
        }

        if (si.hasAmbiguities() && !trans_model.no_ambig_resolve)
        {
            ML_string_distance ml_dist = extendedDistanceInfo.getDistance(i, closestNeig).second;
            si.resolveAmbiguitiesUsingTransitionProbabilities(seqs[closestNeig], ml_dist);
        }
    }

    // UPDATE USING AMBIGUITIES
    if (!trans_model.no_ambiguities)
    {
        for (size_t i = 0; i < numSequences; i++)
        {
            DNA_b128_String &si = seqs[i];
            for (size_t j = i + 1; j < numSequences; j++)
            {
                if (si.hasAmbiguities() || seqs[j].hasAmbiguities())
                {
                    RawDistance sd = extendedDistanceInfo.getDistance(i, j).first;
                    ML_string_distance ml_dist = extendedDistanceInfo.getDistance(i, j).second;
                    sd = correctAmbiguities(sd, ml_dist, si, seqs[j]);

                    ml_dist = computeML(sd);
                    dm.setDistance(i, j, ml_dist.distance);
                }
            }
        }
    }
}

// Shared by fillMatrixRow_{Hamming,JC,K2P,TN93} - see the
// consolidation comment above for why the row-streaming family unifies
// across all four models where the full-matrix family doesn't.
template <class RawDistance, class ComputeRawFn, class ComputeFinalFn>
void fillMatrixRowGeneric(StrFloRow &dm, std::vector<DNA_b128_String> &seqs, size_t row, bool mem_eff_flag,
                          ComputeRawFn computeRaw, ComputeFinalFn computeFinal)
{
    const size_t numSequences = seqs.size();
    dm.resize(numSequences);

    DNA_b128_String &si = seqs[row];
    dm.setDistance(row, 0);

    // Non-memory-efficient (regular row streaming, e.g.
    // --output-format=binary): only the upper triangle (j>row) is
    // computed; j<row is explicitly zeroed, since those distances were
    // already produced when this pair was computed as an earlier row.
    // Memory-efficient (--memory-efficient): every row is
    // self-contained, so the whole range [0,numSequences) is
    // recomputed from scratch, including a harmless redundant
    // self-comparison at j==row (computeRaw(si,si) evaluates to ~0,
    // the same value dm.setDistance(row,0) above already set).
    const size_t startCol = mem_eff_flag ? 0 : row + 1;
    if (!mem_eff_flag)
    {
        for (size_t j = 0; j < row; j++)
        {
            dm.setDistance(j, 0);
        }
    }

    for (size_t j = startCol; j < numSequences; j++)
    {
        RawDistance sd = computeRaw(si, seqs[j]);
        dm.setDistance(j, computeFinal(sd));
    }
}

} // namespace

void
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
fillMatrix_Hamming(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs, sequence_translation_model trans_model)
{

    const size_t numSequences = seqs.size();
    const size_t strlen = seqs[0].getNumChars();

    dm.resize(numSequences);

    // extended information used to compute ambiguities
    static DistanceMatrix<int, simple_string_distance, empty_Data_init<int>, empty_Data_printOn<int>,
                          empty_Data_init<simple_string_distance>, empty_Data_printOn<simple_string_distance>>
        extendedDistanceInfo(numSequences);
    extendedDistanceInfo.resize(numSequences);

    //
    // The loop will resolve the ambiguities according to the translation model
    //
    for (size_t i = 0; i < numSequences; i++)
    {
        dm.setDistance(i, i, 0);
        size_t closestNeig = i;
        float closestDist = FLT_MAX;

        DNA_b128_String &si = seqs[i];
        // only if it has ambiguities does it need to be resolved
        // start by checking the allready computed distances for si
        if (si.hasAmbiguities())
        {
            for (size_t k = 0; k < i; k++)
            {
                float dist = dm.getDistance(k, i);
                if (dist < closestDist && dist >= 0)
                {
                    closestDist = dist;
                    closestNeig = k;
                }
            }
        }

        // compute the remaining distances for si
        for (size_t j = i + 1; j < numSequences; j++)
        {
            // compute distance without using the ambiguities
            simple_string_distance sd = DNA_b128_String::computeDistance(si, seqs[j]);
            float hamdist = compute_Hamming_distance(sd);

            dm.setDistance(i, j, hamdist);
            extendedDistanceInfo.setDistance(i, j, sd);

            if (hamdist < closestDist && hamdist >= 0)
            { // update the closest neighbor
                closestDist = hamdist;
                closestNeig = j;
            }
        }
        // all distances for si have been compted

        // Resolve the ambiguities according to the closest neighbor
        if (si.hasAmbiguities() && !trans_model.no_ambig_resolve)
        {
            si.resolveAmbiguities(seqs[closestNeig]);
        }
    }

    // Update the computed distances with the ambiguities
    if (!trans_model.no_ambiguities)
    {
        for (size_t i = 0; i < numSequences; i++)
        {
            DNA_b128_String &si = seqs[i];
            for (size_t j = i + 1; j < numSequences; j++)
            {
                if (si.hasAmbiguities() || seqs[j].hasAmbiguities())
                {
                    simple_string_distance sd = extendedDistanceInfo.getDistance(i, j);
                    ML_string_distance ml_dist = compute_JC(static_cast<int>(strlen), sd);
                    sd = DNA_b128_String::correctDistanceWithAmbiguitiesUsingTransitionProbabilities(sd, ml_dist, si,
                                                                                                     seqs[j]);

                    // hamming
                    dm.setDistance(i, j, compute_Hamming_distance(sd));
                }
            }
        }
    }
}

void fillMatrix_JC(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs, sequence_translation_model trans_model)
{

    const size_t strlen = seqs[0].getNumChars();

    fillMatrixML<simple_string_distance>(
        dm, seqs, trans_model,
        [](const DNA_b128_String &a, const DNA_b128_String &b) { return DNA_b128_String::computeDistance(a, b); },
        [strlen](simple_string_distance sd) { return compute_JC(static_cast<int>(strlen), sd); },
        [](simple_string_distance sd, ML_string_distance ml, const DNA_b128_String &a, const DNA_b128_String &b) {
            return DNA_b128_String::correctDistanceWithAmbiguitiesUsingTransitionProbabilities(sd, ml, a, b);
        });
}

void fillMatrix_K2P(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs, sequence_translation_model trans_model)
{

    const size_t strlen = seqs[0].getNumChars();

    fillMatrixML<simple_string_distance>(
        dm, seqs, trans_model,
        [](const DNA_b128_String &a, const DNA_b128_String &b) { return DNA_b128_String::computeDistance(a, b); },
        [strlen, &trans_model](simple_string_distance sd) {
            if (trans_model.no_tstvratio)
            {
                return compute_K2P(static_cast<int>(strlen), sd);
            }
            return compute_K2P_fixratio(static_cast<int>(strlen), sd, trans_model.tstvratio);
        },
        [](simple_string_distance sd, ML_string_distance ml, const DNA_b128_String &a, const DNA_b128_String &b) {
            return DNA_b128_String::correctDistanceWithAmbiguitiesUsingTransitionProbabilities(sd, ml, a, b);
        });
}

void fillMatrix_TN93(StrDblMatrix &dm, std::vector<DNA_b128_String> &seqs, DNA_b128_String::base_frequences freqs,
                     sequence_translation_model trans_model)
{

    const size_t strlen = seqs[0].getNumChars();

    fillMatrixML<TN_string_distance>(
        dm, seqs, trans_model,
        [](const DNA_b128_String &a, const DNA_b128_String &b) {
            return DNA_b128_String::computeTAMURANEIDistance(a, b);
        },
        [strlen, freqs, &trans_model](TN_string_distance tn) {
            // compute_Tamura_Nei{,_fixratio}() must keep int strlen/numAs/
            // numCs/numGs/numTs (lint_plan.md's Phase 2 investigation:
            // strlen participates in a signed subtraction internally) -
            // cast once here rather than changing their signatures.
            const int strlenAsInt = static_cast<int>(strlen);
            const int numAs = static_cast<int>(freqs.num_As_);
            const int numCs = static_cast<int>(freqs.num_Cs_);
            const int numGs = static_cast<int>(freqs.num_Gs_);
            const int numTs = static_cast<int>(freqs.num_Ts_);
            if (trans_model.no_tstvratio)
            {
                return compute_Tamura_Nei(strlenAsInt, tn, numAs, numCs, numGs, numTs);
            }
            return compute_Tamura_Nei_fixratio(strlenAsInt, tn, numAs, numCs, numGs, numTs, trans_model.tstvratio,
                                               trans_model.pyrtvratio);
        },
        [&trans_model](TN_string_distance tn, ML_string_distance ml, const DNA_b128_String &a,
                       const DNA_b128_String &b) {
            if (trans_model.no_transition_probs)
            {
                return DNA_b128_String::correctDistanceWithAmbiguitiesUsingBackgroundFrequences(tn, a, b);
            }
            return DNA_b128_String::correctDistanceWithAmbiguitiesUsingTransitionProbabilities(tn, ml, a, b);
        });
}

//-------- mehmood warka dang---------------

void fillMatrixRow(StrFloRow &dm, std::vector<DNA_b128_String> &seqs, sequence_translation_model trans_model,
                   size_t row, bool mem_eff_flag)
{

    const size_t numSequences = seqs.size();
    dm.resize(numSequences);

    // base frequences
    DNA_b128_String::base_frequences freqs = seqs[0].getBaseFrequences();
    for (size_t i = 1; i < numSequences; i++)
    {
        DNA_b128_String::base_frequences tmpfreqs = seqs[i].getBaseFrequences();
        freqs.num_As_ += tmpfreqs.num_As_;
        freqs.num_Cs_ += tmpfreqs.num_Cs_;
        freqs.num_Gs_ += tmpfreqs.num_Gs_;
        freqs.num_Ts_ += tmpfreqs.num_Ts_;
        freqs.num_unknowns_ += tmpfreqs.num_unknowns_;
        freqs.num_ambiguities_ += tmpfreqs.num_ambiguities_;
    }

    // compute ambiguities probabilities
    if (!trans_model.no_ambiguities)
    {
        if (trans_model.use_base_freqs)
        {
            for (size_t i = 0; i < numSequences; i++)
            {
                seqs[i].calcAmbiguityProbabilities(freqs.num_As_, freqs.num_Cs_, freqs.num_Gs_, freqs.num_Ts_);
            }
        }
        else
        {
            for (size_t i = 0; i < numSequences; i++)
            {
                seqs[i].calcAmbiguityProbabilitiesUNIFORM();
            }
        }
    }

    // CALL DISTANCE COMPUTATION
    switch (trans_model.model)
    {
    case HAMMING_DISTANCE:
        fillMatrixRow_Hamming(dm, seqs, trans_model, row, mem_eff_flag);
        break;
    case JC:
        fillMatrixRow_JC(dm, seqs, trans_model, row, mem_eff_flag);
        break;
    case K2P:
        fillMatrixRow_K2P(dm, seqs, trans_model, row, mem_eff_flag);
        break;
    case TN93:
        fillMatrixRow_TN93(dm, seqs, freqs, trans_model, row, mem_eff_flag);
        break;
    default:
        THROW_EXCEPTION("Non handled model: " << trans_model.model);
    }
}

void fillMatrixRow_Hamming(StrFloRow &dm, std::vector<DNA_b128_String> &seqs, sequence_translation_model trans_model,
                           size_t row, bool mem_eff_flag)
{
    fillMatrixRowGeneric<simple_string_distance>(
        dm, seqs, row, mem_eff_flag,
        [](const DNA_b128_String &a, const DNA_b128_String &b) { return DNA_b128_String::computeDistance(a, b); },
        [](simple_string_distance sd) { return compute_Hamming_distance(sd); });
}

void fillMatrixRow_JC(StrFloRow &dm, std::vector<DNA_b128_String> &seqs, sequence_translation_model trans_model,
                      size_t row, bool mem_eff_flag)
{

    const size_t strlen = seqs[0].getNumChars();

    fillMatrixRowGeneric<simple_string_distance>(
        dm, seqs, row, mem_eff_flag,
        [](const DNA_b128_String &a, const DNA_b128_String &b) { return DNA_b128_String::computeDistance(a, b); },
        [strlen](simple_string_distance sd) { return compute_JC(static_cast<int>(strlen), sd).distance; });
}

void fillMatrixRow_K2P(StrFloRow &dm, std::vector<DNA_b128_String> &seqs, sequence_translation_model trans_model,
                       size_t row, bool mem_eff_flag)
{

    const size_t strlen = seqs[0].getNumChars();

    fillMatrixRowGeneric<simple_string_distance>(
        dm, seqs, row, mem_eff_flag,
        [](const DNA_b128_String &a, const DNA_b128_String &b) { return DNA_b128_String::computeDistance(a, b); },
        [strlen, &trans_model](simple_string_distance sd) {
            if (trans_model.no_tstvratio)
            {
                return compute_K2P(static_cast<int>(strlen), sd).distance;
            }
            return compute_K2P_fixratio(static_cast<int>(strlen), sd, trans_model.tstvratio).distance;
        });
}

void fillMatrixRow_TN93(StrFloRow &dm, std::vector<DNA_b128_String> &seqs, DNA_b128_String::base_frequences freqs,
                        sequence_translation_model trans_model, size_t row, bool mem_eff_flag)
{

    const size_t strlen = seqs[0].getNumChars();

    fillMatrixRowGeneric<TN_string_distance>(
        dm, seqs, row, mem_eff_flag,
        [](const DNA_b128_String &a, const DNA_b128_String &b) {
            return DNA_b128_String::computeTAMURANEIDistance(a, b);
        },
        [strlen, freqs, &trans_model](TN_string_distance tn) {
            // Same must-stay-int rationale as fillMatrix_TN93 above.
            const int strlenAsInt = static_cast<int>(strlen);
            const int numAs = static_cast<int>(freqs.num_As_);
            const int numCs = static_cast<int>(freqs.num_Cs_);
            const int numGs = static_cast<int>(freqs.num_Gs_);
            const int numTs = static_cast<int>(freqs.num_Ts_);
            if (trans_model.no_tstvratio)
            {
                return compute_Tamura_Nei(strlenAsInt, tn, numAs, numCs, numGs, numTs).distance;
            }
            return compute_Tamura_Nei_fixratio(strlenAsInt, tn, numAs, numCs, numGs, numTs, trans_model.tstvratio,
                                               trans_model.pyrtvratio)
                .distance;
        });
}

//----------- khatam sho janana--------------
