#include "fastphylo/io/BinaryDmOutputStream.hpp"
#include "fastphylo/core/file_utils.hpp"
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

BinaryDmOutputStream::BinaryDmOutputStream(char *filename) : DataOutputStream(filename) {
  if (filename != nullptr) {
    writeToCout = false;
    // The base DataOutputStream(filename) constructor already opened fp
    // via fopen() (open_write_file()); we don't need it (we write
    // through ofs instead), so close it properly.
    fclose(fp);
    fp = nullptr;
    file_was_opened = false;
    file_stream.reset(open_write_binary(filename));
    ofs = file_stream.get();
  } else {
    ofs = &cout;
    writeToCout = true;
  }
}

void BinaryDmOutputStream::print(StrDblMatrix &dm) {
  // One write() per row instead of one per float - same batching
  // rationale as PhylipDmOutputStream.cpp's printPHYLIPfastSD().
  // std::vector<float> is contiguous/tightly-packed, so this writes the
  // exact same bytes in the exact same order as before, just in
  // numNodes calls instead of numNodes*(numNodes+1)/2.
  const size_t numNodes = dm.getSize();
  vector<float> row;
  row.reserve(numNodes);
  for (size_t i = 0; i < numNodes; i++) {
    row.clear();
    for (size_t j = i; j < numNodes; j++) {
      float f = dm.getDistance(i, j);
      if (!isfinite(f))
        f = -1.0;
      row.push_back(f);
    }
    ofs->write(reinterpret_cast<char *>(row.data()), row.size() * sizeof(float));
  }
}

// fastdist's bootstrap-streaming path: writes one row at a time as it's
// computed, rather than building the whole matrix first.
void BinaryDmOutputStream::printRow(StrFloRow &dm, string name, int row, bool mem_eff_flag) {
  int entriesPerRow = dm.getColumns();

  for (size_t j = row; j < static_cast<size_t>(entriesPerRow); j++) {
    float f = dm.getDistance(j);

    if (!isfinite(f)) {
      f = -1.0;
    }

    ofs->write(reinterpret_cast<char *>(&f), sizeof f);
  }
}

void BinaryDmOutputStream::printHeader(size_t numNodes) {
  // converter variable is needed for running the binary output/input
  // also on 64-bit systems
  string tag = "FASTPHYLO 1";
  ofs->write(tag.c_str(), tag.length());
  long converter = numNodes;
  ofs->write(reinterpret_cast<char *>(&converter), sizeof(converter));
  for (size_t i = 0; i < numNodes; ++i) {
    string name = m_names[i];
    ofs->write(name.c_str(), name.length());
    ofs->write(":", 1);
  }
}

void BinaryDmOutputStream::printStartRun(std::vector<std::string> &names, std::string &runId, Extrainfos &extrainfos) {
  m_names = names;
}
