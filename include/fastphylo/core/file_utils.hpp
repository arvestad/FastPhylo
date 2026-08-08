//--------------------------------------------------
//
// File: file_utils.hpp
//
//--------------------------------------------------

#pragma once

#include <cstdio>
#include <fstream>
#include <istream>

// legacy_error_handling_plan.md Finding 3/Phase 4: this used to declare
// 14 functions, mostly char*-first with a std::string overload as a
// thin .c_str()-forwarding afterthought. A direct audit (redone after
// Phase 0 deleted several files that were file_utils' only remaining
// callers of open_write_stream/open_read_stream) found only these 3
// with any live caller left anywhere in the tree; everything else -
// the "interactive" prompt-and-retry family, open_write_stream,
// open_read_stream, file_exists (only ever called by the two deleted
// "interactive" functions), skipWhiteSpace(FILE*), skipUntil,
// appendToken, appendUntil - was dead and has been removed.
FILE *open_write_file(const char *fname);

// Added by Mehmood Khan Malagori; email: malagori@kth.se
std::ofstream *open_write_binary(const char *fname);

void skipWhiteSpace(std::istream &in);
