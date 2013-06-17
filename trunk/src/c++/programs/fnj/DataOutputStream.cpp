#include "DataOutputStream.hpp"
#include <cstdio>

using namespace std;

DataOutputStream::DataOutputStream(char * filename) {
  file_was_opened = false;
  if (filename==NULL)
    fp = &cout;
  else {
    fout.open(filename, ofstream::out);
    if (!fout.good()) {
      fout.close();
      fout.clear();
      THROW_EXCEPTION("File doesn't exist: \"" << filename << "\"");
      }
    file_was_opened = true;
    fp = & fout;
  }
}

