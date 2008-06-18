#include "XmlOutputStream.hpp"
#include <cstdio>
#include <libxml/xmlreader.h>

using namespace std;

XmlOutputStream::XmlOutputStream(char * filename = 0 ) : DataOutputStream(filename) 
{
  LIBXML_TEST_VERSION
   *fp << "<?xml version=\"1.0\"?>" << std::endl << " <root>" << std::endl << "  <runs>" <<  std::endl ;
};

XmlOutputStream::~XmlOutputStream() 
{
  *fp << "  </runs>" << std::endl << " </root>" <<  std::endl ;
};

void
XmlOutputStream::print( tree2int_map & tree2count, bool noCounts ) 
{

  tree2int_map::iterator iter = tree2count.begin();
  for( ; iter!=tree2count.end() ; ++iter){
    *fp << "   <run>" <<  std::endl
        << "    <count>"  << (*iter).second 
        << "</count>"  <<  std::endl 
        << "    <tree>" << (*iter).first
        << "    </tree>" << std::endl
        << "   </run>"  << std::endl;
  }
}

