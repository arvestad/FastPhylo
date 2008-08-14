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
XmlOutputStream::print( tree2int_map & tree2count, bool noCounts, std::string & runId, std::vector<std::string> & names, Extrainfos & extrainfos ) 
{
  Extrainfos::iterator it;
  std::vector<std::string>::iterator it2;
  *fp << "   <run id=\""   <<  runId <<      "\" dim=\"" <<  names.size()  << "\">" <<  std::endl 
      << "    <identities>" <<  std::endl;
  it=extrainfos.begin();
  for (it2=names.begin() ; it2 != names.end()  ;  ++it2 )
    {  
      if ( it != extrainfos.end() && it->size() > 0 ) { 
	*fp << "     <identity name=\""  <<  *it2  <<  "\">" << *it << "     </identity>" <<  std::endl;         
        ++it;
      } 
      else {
	*fp << "     <identity name=\""  <<  *it2  <<  "\"/>" <<  std::endl;         
      }
    }

  *fp  << "    </identities>" <<  std::endl;

  tree2int_map::iterator iter = tree2count.begin();
  for( ; iter!=tree2count.end() ; ++iter){
    *fp  << "    <tree>" <<  std::endl
<< "     <count>"  << (*iter).second 
        << "</count>"  <<  std::endl 
        << "     <newick-xml>" << (*iter).first
	 << "</newick-xml>" << std::endl
         << "     <newick>";


    ostringstream oss;
    oss <<  (*iter).first;
    printNewick( fp , oss.str() );
    *fp   << ";</newick>"  << std::endl
    << "    </tree>"  << std::endl;
    }
   *fp     << "   </run>"  << std::endl;
}

void XmlOutputStream::printNewick(std::ostream * fp , std::string s ) {
  xmlDocPtr doc;
  xmlNode *root = NULL;
  doc = xmlReadMemory(s.c_str(), s.size(), "", NULL, 0);
  if (doc == NULL) {
    std::cerr << "internal parse error" << std::endl; 
    return;
  }
  root = xmlDocGetRootElement(doc);
  printNewickNode(fp , root);
  xmlFreeDoc(doc);
 
  return;
}

void XmlOutputStream::printNewickNode(std::ostream * fp , xmlNode * node)
{
  xmlNode *n = NULL;
 
   bool firsttime = true;
   for (n = node; n; n = n->next) {
     if (n->type == XML_ELEMENT_NODE) {
       if ( firsttime ) {
	 firsttime = false;
       } else {
	 *fp << ",";
       }
       if ( xmlStrEqual (n->name, (const xmlChar *)"branch" ) ) {
         *fp << "("; 
	 printNewickNode( fp, n->children );
         *fp << ")"; 
       } 
       if ( xmlStrEqual (n->name, (const xmlChar *)"leaf" ) ) {
	 *fp << xmlNodeGetContent(n);
       } 
     }
   }
   return;
}

