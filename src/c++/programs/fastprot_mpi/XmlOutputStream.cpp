#include "XmlOutputStream.hpp"
#include <cstdio>
#include <string>


using namespace std;

XmlOutputStream::XmlOutputStream(char * filename = 0 ) : DataOutputStream(filename) 
{
  fprintf(fp,"<?xml version=\"1.0\"?>\n<root>\n <runs>\n");
};

XmlOutputStream::~XmlOutputStream() 
{
  fprintf(fp," </runs>\n</root>\n");
};

void
XmlOutputStream::print( StrDblMatrix & dm ) 
{
  //printPHYLIPfast(dm, fp , true );
  printPHYLIPfastSD(dm, fp , true, false );
}

void
XmlOutputStream::printSD( StrDblMatrix & dm ) 
{
  printPHYLIPfastSD(dm, fp , false, true );
}

void
XmlOutputStream::printStartRun( std::vector<string> & names, std::string & runId, Extrainfos &extrainfos ) 
{
  fprintf(fp,"  <run id=\"%s\" dim=\"%i\">\n   <identities>\n", runId.c_str(), names.size() );
      for(size_t namei=0 ; namei<names.size() ; namei++ )
	{
	  // This only works if we constrain the input by a schema to not have "<", "&" and such.
          // Otherwise we need to use xmlEncodeSpecialChars(xmlDocPtr doc,  const xmlChar * input)
	  //	  xmlChar * str = xmlEncodeSpecialChars( 0 ,( const xmlChar * ) names[namei].c_str()  );
	  fprintf(fp,"    <identity name=\"%s\"", names[namei].c_str() );

	  if ( extrainfos.size() > namei  && extrainfos[namei].size() > 0 ) {
	    fprintf(fp,">\n     %s\n    </identity>\n", (char *) extrainfos[namei].c_str() );
          } 
          else {
	    fprintf(fp,"/>\n");
         } 

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



