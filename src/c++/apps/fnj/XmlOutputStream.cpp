#include "XmlOutputStream.hpp"
#include <cstdio>

using namespace std;

XmlOutputStream::XmlOutputStream(char * filename) : DataOutputStream(filename) {
  *fp << "<?xml version=\"1.0\"?>" << std::endl << " <root>" << std::endl << "  <runs>" <<  std::endl;
}

XmlOutputStream::~XmlOutputStream() {
  *fp << "  </runs>" << std::endl << " </root>" <<  std::endl;
}

void XmlOutputStream::print( tree2int_map & tree2count, bool printCounts,string & runId,vector<string> & names, Extrainfos & extrainfos) {
  Extrainfos::iterator it;
  *fp << "   <run id=\""   <<  runId <<      "\" dim=\"" <<  names.size()  << "\">" <<  std::endl
      << "    <identities>" <<  std::endl;
  it=extrainfos.begin();
  for (const string &name : names)
    if ( it != extrainfos.end() && it->size() > 0 ) {
      *fp << "     <identity name=\""  <<  name  <<  "\">" << *it << "     </identity>" <<  std::endl;
      ++it;
    }
    else
      *fp << "     <identity name=\""  <<  name  <<  "\"/>" <<  std::endl;
  *fp  << "    </identities>" <<  std::endl;
  for (const auto &entry : tree2count) {
    ostringstream oss;
    *fp  << "    <tree>" <<  std::endl
      << "     <count>"  << entry.second
      << "</count>"  <<  std::endl;
    bool oldVal = xmlPrint;
    xmlPrint = true;
    *fp  << "     <newick-xml>" <<  entry.first <<  "</newick-xml>" << std::endl;
    xmlPrint = false;
    *fp  << "     <newick>" <<  entry.first <<  "</newick>" << std::endl
      << "    </tree>"  << std::endl;
    xmlPrint = oldVal;
  }
  *fp     << "   </run>"  << std::endl;
}
