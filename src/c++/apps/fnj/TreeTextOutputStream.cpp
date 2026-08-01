#include "TreeTextOutputStream.hpp"
#include <cstdio>

using namespace std;

void TreeTextOutputStream::print(tree2int_map & tree2count, bool printCounts,string & runId,vector<string> & names, Extrainfos & extrainfos) {
  //OUTPUT THE TREES
  for (const auto &entry : tree2count) {
    if(printCounts)
      *fp << entry.second << "  " << entry.first << endl;
    else
      *fp << entry.first << endl;
  }
}

