//--------------------------------------------------
//
// File: Exception.hpp
//
//--------------------------------------------------
#pragma once

#include <sstream>
#include <stdexcept>

// legacy_error_handling_plan.md Finding 2/Phase 3: this used to define a
// hand-rolled Exception class (file/function/line fields, its own
// printOn()-based formatting, double-inheriting a Java-style Object base
// and std::exception without overriding what()). A direct audit found
// nothing anywhere in the tree ever read those fields, caught Exception
// structurally, or used its dead addToStackTrace()/CATCH_RETHROW()
// mechanism - every one of THROW_EXCEPTION's 57 call sites just streams
// a message. std::runtime_error already does exactly that, correctly,
// with no custom class needed - THROW_EXCEPTION now builds the same
// file/line/function-prefixed message text and throws that directly.
#define THROW_EXCEPTION(MES) \
  { std::ostringstream out; out << __FILE__ << ":" << __LINE__ << " (" << __FUNCTION__ << "): " << MES; \
    throw std::runtime_error(out.str()); }
