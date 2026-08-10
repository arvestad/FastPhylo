#include "fastphylo/protein/ProtSeqCompare.hpp"
#include "fastphylo/protein/ProtSeqCode.hpp"
#include <algorithm>

// Real SSE2 on x86; simde's NEON-backed translation everywhere else
// (e.g. Apple Silicon) - same pattern as DNA_b128/sse2_wrapper.h.
// SIMDE_ENABLE_NATIVE_ALIASES makes simde define the plain _mm_* names
// this file calls, instead of simde_mm_*.
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
#include <emmintrin.h>
#else
#define SIMDE_ENABLE_NATIVE_ALIASES
#include <simde/x86/sse2.h>
#endif

// __builtin_popcount is a GCC/Clang compiler builtin (works identically
// on the ARM/simde path too - it isn't SSE-specific, every GCC/Clang
// target supports it), not part of any standard header, and MSVC does
// not provide it at all. __popcnt (<intrin.h>) is MSVC's equivalent
// intrinsic. clang-cl defines _MSC_VER but still supports the GCC
// builtin natively, so it is excluded from the MSVC branch.
#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>
#endif

namespace ProtSeqCode
{

namespace
{

inline unsigned int popcount16(unsigned int mask)
{
#if defined(_MSC_VER) && !defined(__clang__)
    return __popcnt(mask);
#else
    return static_cast<unsigned int>(__builtin_popcount(mask));
#endif
}

// Count mismatching bytes between a and b over exactly n bytes.
// 16 bytes/iteration via SSE2 pcmpeqb+pmovmskb+popcount (NEON on
// arm64 via simde translation); scalar loop for the remainder.
std::size_t count_mismatches_exact(const std::uint8_t *a, const std::uint8_t *b, std::size_t n)
{
    std::size_t mismatches = 0;
    std::size_t i = 0;

    for (; i + 16 <= n; i += 16)
    {
        __m128i va = _mm_loadu_si128(reinterpret_cast<const __m128i *>(a + i));
        __m128i vb = _mm_loadu_si128(reinterpret_cast<const __m128i *>(b + i));
        __m128i eq = _mm_cmpeq_epi8(va, vb);
        auto eq_mask = static_cast<unsigned int>(_mm_movemask_epi8(eq));
        unsigned int mismatch_mask = (~eq_mask) & 0xFFFFU;
        mismatches += static_cast<std::size_t>(popcount16(mismatch_mask));
    }

    for (; i < n; i++)
    {
        if (a[i] != b[i])
        {
            mismatches++;
        }
    }

    return mismatches;
}

} // anonymous namespace

std::size_t count_mismatches(const std::uint8_t *s1, std::size_t len1, const std::uint8_t *s2, std::size_t len2)
{
    std::size_t common = std::min(len1, len2);
    std::size_t mismatches = count_mismatches_exact(s1, s2, common);
    if (len1 > common)
    {
        mismatches += (len1 - common);
    }
    return mismatches;
}

double count_id_fraction(const std::uint8_t *s1, std::size_t len1, const std::uint8_t *s2, std::size_t len2)
{
    std::size_t mismatches = count_mismatches(s1, len1, s2, len2);
    double matches = static_cast<double>(len1) - static_cast<double>(mismatches);
    return matches / static_cast<double>(len1);
}

std::vector<std::size_t> count_replacement_tally(const std::uint8_t *s1, std::size_t len1, const std::uint8_t *s2,
                                                 std::size_t len2)
{
    std::vector<std::size_t> counts(NUM_CANONICAL_AA * NUM_CANONICAL_AA, 0);
    std::size_t n = std::min(len1, len2);
    for (std::size_t i = 0; i < n; i++)
    {
        std::uint8_t c1 = s1[i];
        std::uint8_t c2 = s2[i];
        if (is_canonical_aa(c1) && is_canonical_aa(c2))
        {
            counts[(c1 * NUM_CANONICAL_AA) + c2]++;
        }
    }
    return counts;
}

} // namespace ProtSeqCode
