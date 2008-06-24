#include "XmlOutputStream.hpp"
#include <cstdio>
#include <string>
#include <libxml/xmlreader.h>

using namespace std;

XmlOutputStream::XmlOutputStream(char * filename = 0 ) : DataOutputStream(filename) 
{
  LIBXML_TEST_VERSION
  fprintf(fp,"<?xml version=\"1.0\"?>\n<root>\n <runs>\n");
};

XmlOutputStream::~XmlOutputStream() 
{
  fprintf(fp," </runs>\n</root>\n");
};

void
XmlOutputStream::print( StrDblMatrix & dm ) 
{
  printPHYLIPfast(dm, fp , true );
}


void
XmlOutputStream::printStartRun( std::vector<string> & names ) 
{
  fprintf(fp,"  <run dim=\"%i\">\n   <identities>\n", names.size() );
      for(size_t namei=0 ; namei<names.size() ; namei++ )
	{
	  // This only works if we constrain the input by a schema to not have "<", "&" and such.
          // Otherwise we need to use xmlEncodeSpecialChars(xmlDocPtr doc,  const xmlChar * input)
	  //	  xmlChar * str = xmlEncodeSpecialChars( 0 ,( const xmlChar * ) names[namei].c_str()  );
	          fprintf(fp,"    <identity name=\"%s\"/>\n", names[namei].c_str() );
		  // fprintf(fp,"     <entry>%s</entry>\n", ( char * ) str );
		  //  xmlFree(str);
	}
  fprintf(fp,"   </identities>\n   <dms>\n");
}


void
XmlOutputStream::printEndRun() 
{
  fprintf(fp,"   </dms>\n  </run>\n");
}



