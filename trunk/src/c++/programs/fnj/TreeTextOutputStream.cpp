#include "TreeTextOutputStream.hpp"
#include <cstdio>

using namespace std;

void TreeTextOutputStream::print(tree2int_map & tree2count, bool printCounts,string & runId,vector<string> & names, Extrainfos & extrainfos) {
  //OUTPUT THE TREES
  tree2int_map::iterator iter = tree2count.begin();
  for(; iter!=tree2count.end(); iter++) {
    if(printCounts)
      *fp << iter->second << "  " << iter->first << endl;
    else
      *fp << iter->first << endl;
  }
}

