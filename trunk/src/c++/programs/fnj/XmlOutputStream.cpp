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

    ostringstream oss;
    oss <<  (*iter).first;

    *fp  << "    <tree>" <<  std::endl
<< "     <count>"  << (*iter).second 
        << "</count>"  <<  std::endl 
	 << "     <newick-xml>";

    printNewick( fp , oss.str(), true );
  *fp 	 << "</newick-xml>" << std::endl
         << "     <newick>";

    printNewick( fp , oss.str(), false );
    *fp   << ";</newick>"  << std::endl
    << "    </tree>"  << std::endl;
    }
   *fp     << "   </run>"  << std::endl;
}

void XmlOutputStream::printNewick(std::ostream * fp , std::string s, bool printXml ) {
  xmlDocPtr doc;
  xmlNode *root = NULL;

  std:string rootstr = "<root>" + s + "</root>"; 
  doc = xmlReadMemory(rootstr.c_str(), rootstr.size(), "", NULL, 0);
  if (doc == NULL) {
    std::cerr << "internal parse error" << std::endl; 
    return;
  }
  root = xmlDocGetRootElement(doc);
  printNewickNode(fp , root, printXml);
  xmlFreeDoc(doc); 
  return;
}

void XmlOutputStream::printNewickNode(std::ostream * fp , xmlNode * node, bool printXml)
{
  xmlNode *n = NULL;
 
   bool firsttime = true;

   for (n = node; n; n = n->next) {
     if (n->type == XML_ELEMENT_NODE  ) {
       //       *fp << "n->name=" << (char * ) n->name << std::endl;
       if (  printXml == false && ( xmlStrEqual (n->name, (const xmlChar *)"branch") ||  xmlStrEqual (n->name, (const xmlChar *)"leaf")) ) {
         if ( firsttime == true ) {
            firsttime = false;
         } else {
	     *fp << ",";
         }
       }
       if (  printXml == true && xmlStrEqual (n->name, (const xmlChar *)"branch")  ) {
	   *fp << "<branch";
             if ( n->next && n->next->type == XML_ELEMENT_NODE && xmlStrEqual (n->next->name, (const xmlChar *)"length")  ) {
	       *fp << " length=\"";
	       xmlChar * leng = xmlNodeGetContent(n->next);
	       *fp <<  leng << "\"" ;
	       xmlFree(leng);
	     }
   	     *fp << ">";
	     printNewickNode( fp, n->children, printXml );
	   *fp << "</branch>";
       }
       if (  printXml == true && xmlStrEqual (n->name, (const xmlChar *)"leaf")  ) {
	   *fp << "<leaf";
	   if (  n->children && n->children->next && n->children->next->type == XML_ELEMENT_NODE &&
  xmlStrEqual (n->children->next->name, (const xmlChar *)"length")   ) {
	       *fp << " length=\"";
	       xmlChar * leng = xmlNodeGetContent(n->children->next);
	       *fp <<  leng << "\"" ;
	       xmlFree(leng);
	     }
   	     *fp << ">";
	   if (  n->children && n->children->type == XML_TEXT_NODE  ) {
	     xmlChar * text = xmlNodeGetContent(n->children);
	     *fp <<  text  ;
	     xmlFree(text);
	   }
	   *fp << "</leaf>";
      }
      if ( printXml == false && xmlStrEqual (n->name, (const xmlChar *)"length" ) ) {
           xmlChar * len = xmlNodeGetContent(n);
	   *fp <<  ":" << len ;
           xmlFree(len);
       }
       if ( printXml == false && xmlStrEqual (n->name, (const xmlChar *)"branch" ) ) {
         *fp << "("; 
	 printNewickNode( fp, n->children, printXml );
         *fp << ")"; 
       } 
       if ( xmlStrEqual (n->name, (const xmlChar *)"root" ) ) {
	 printNewickNode( fp, n->children, printXml );
       } 
       if ( printXml == false && xmlStrEqual (n->name, (const xmlChar *)"leaf" ) && n->children != NULL ) {
         xmlChar * species = xmlNodeGetContent(n->children);
       	 *fp << species;
         xmlFree(species);
         if ( n->children->next != NULL  && 
           n->children->next->type == XML_ELEMENT_NODE &&
           xmlStrEqual (n->children->next->name, (const xmlChar *)"length" ) ) {
           xmlChar * length = xmlNodeGetContent(n->children->next);
	   *fp <<  ":" << length;
           xmlFree(length);
	 }
       } 
     }
   }
   return;
}

