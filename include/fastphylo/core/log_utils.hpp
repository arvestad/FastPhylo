//--------------------------------------------------
//
// File: log_utils.hpp
//
//--------------------------------------------------
#pragma once

#include <iostream>
#include <cassert>

// legacy_error_handling_plan.md Finding 5 / modern_cpp_core_cleanup_plan.md:
// this used to define 11 macros (PRINT/PRINT_V/PRINT_TIME/PRINT_EXP/LINE/
// SEPARATOR - all dead printf-debugging leftovers with zero real call
// sites; MEM_CHECK/PROG_ERROR/USER_ERROR/USER_WARNING - a second,
// uncoordinated error-reporting system that called exit() directly,
// bypassing THROW_EXCEPTION's throw/unwind/catch entirely). All of
// those are gone now: USER_ERROR/PROG_ERROR call sites became
// THROW_EXCEPTION (identical streaming syntax, but unwinds properly -
// exit() skips destructors of in-scope objects, so a mid-write
// std::ofstream could lose buffered data); USER_WARNING call sites
// became direct std::cerr writes (no control-flow behavior needed
// wrapping); MEM_CHECK's 2 sites became a direct null-check + abort()
// (still not an exception - minimizing further allocation while
// already handling OOM is the one part of the old design worth
// keeping). Only ASSERT_EQ remains, still genuinely used.
#ifndef ASSERT_EQ
#ifndef NDEBUG
#define ASSERT_EQ(EXP1, EXP2)                                                                                          \
    if (!((EXP1) == (EXP2)))                                                                                           \
    {                                                                                                                  \
        std::cerr << "************\nASSERT_EQ FAILED\nfile: " << __FILE__ << "\nfunc: " << __FUNCTION__                \
                  << "\nline: " << __LINE__ << "\n"                                                                    \
                  << #EXP1 << " != " << #EXP2 << "\n"                                                                  \
                  << #EXP1 << " = " << (EXP1) << "\n"                                                                  \
                  << #EXP2 << " = " << (EXP2) << "\n**********" << std::endl;                                          \
        assert((EXP1) == (EXP2));                                                                                      \
    }
#else
#define ASSERT_EQ(EXP1, EXP2)
#endif
#endif
