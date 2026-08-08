//--------------------------------------------------
//
// File: file_utils.cpp
//
//--------------------------------------------------

#include "fastphylo/core/file_utils.hpp"
#include "fastphylo/core/Exception.hpp"

using namespace std;

FILE *
open_write_file(const char *fname) {
  // Binary mode, not "w" (text mode) - on Windows the CRT translates
  // every outgoing '\n' to "\r\n" in text mode, silently corrupting
  // output that is meant to be byte-identical across platforms (this
  // project's own code always writes '\n' explicitly, so it never
  // relied on that translation for correctness).
  FILE *ftmp = fopen(fname, "wb");
  if (ftmp == nullptr) {
    THROW_EXCEPTION("Couldn't open file \"" << fname << "\"");
  }
  return ftmp;
}

// Added by Mehmood Khan Malagori; email: malagori@kth.se
std::ofstream *
open_write_binary(const char *fname) {
  auto *ofs = new ofstream(fname, ios::binary);
  if (!ofs->good()) {
    ofs->close();
    ofs->clear();
    THROW_EXCEPTION("Can't open file " << fname);
  }
  return ofs;
}

void
skipWhiteSpace(std::istream &in) {
  in >> ws;
}
