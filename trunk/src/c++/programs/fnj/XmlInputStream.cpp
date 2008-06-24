#include "XmlInputStream.hpp"
#include <cstdio>

using namespace std;

XmlInputStream::~XmlInputStream()
{ 
  if ( reader ) {
    xmlFreeTextReader(reader);
  }
  xmlCleanupParser();
}

XmlInputStream::XmlInputStream(char * filename = 0)  
{
  //
  //  Re: [xml] Why does "-" read from stdin?
  //  http://mail.gnome.org/archives/xml/2007-February/msg00005.html
  //
  // Comment: The special treatment of the filename "-" in the libxml api, is not  
  // a good designed api, but now when it is there let us use it.

  if ( filename == "-" ) {
    THROW_EXCEPTION("file name \"-\" is not allowed. See \n http://mail.gnome.org/archives/xml/2007-February/msg00005.html \n  ");
  }

  if ( filename == 0 ) {
    filename = ( char * ) "-";
  }
  LIBXML_TEST_VERSION

  reader = xmlReaderForFile(filename,0, XML_PARSE_COMPACT | XML_PARSE_NONET );
  if ( reader == 0 ) { THROW_EXCEPTION("Could not open file"); };

  l.in_root =  false;
  l.in_runs =  false;
  l.in_run = false;
  l.in_identities = false;
  l.row_nr = -1;
  l.entry_nr = -1;

  dmSize = 0;

  xmlRelaxNGParserCtxtPtr parserctxt;
  size_t len = strlen(relaxngstr);
  parserctxt = xmlRelaxNGNewMemParserCtxt(relaxngstr,len);
  xmlRelaxNGSetParserErrors(parserctxt,(xmlRelaxNGValidityErrorFunc) fprintf, (xmlRelaxNGValidityWarningFunc) fprintf, stderr);
  xmlRelaxNGPtr schema = NULL;
  schema = xmlRelaxNGParse(parserctxt);
  xmlRelaxNGFreeParserCtxt(parserctxt);

  if ( xmlTextReaderRelaxNGSetSchema( reader,  schema ) != 0 )  { 
    THROW_EXCEPTION("failed to set relax ng schema");
    exit(EXIT_FAILURE);
  }
}

readstatus  
XmlInputStream::readDM( StrDblMatrix & dm ) {
{ 
    const xmlChar *name, *value;

    bool run_read = false;
    int ret;

    while ( ( ret = xmlTextReaderRead(reader)) == 1 )
      {
 	if ( xmlTextReaderIsValid( reader ) != 1 ) { 
          THROW_EXCEPTION("xml input does not validate");
          exit(EXIT_FAILURE);
        } 

	int depth = xmlTextReaderDepth(reader);
	int type = xmlTextReaderNodeType(reader);
	name = xmlTextReaderConstName(reader);


      	if ( l.in_root && l.in_runs && l.in_run && l.in_dms && l.in_dm &&  l.in_row && 
               depth == 6 &&  xmlStrEqual (name, (const xmlChar *)"entry" )  && type == XML_READER_TYPE_ELEMENT  )
	  {
                l.entry_nr++;
  	        xmlChar * distanceStr = xmlTextReaderReadString(reader);
                float distance =  atof( ( char * ) distanceStr );
                dm.setDistance(l.row_nr,l.entry_nr, distance );  
                xmlFree(distanceStr);
           }

      	if ( l.in_root && l.in_runs && l.in_run && l.in_dms && l.in_dm &&  depth == 5 &&  xmlStrEqual (name, (const xmlChar *)"row" ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_row = true; l.row_nr++; l.entry_nr = -1; continue; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_row = false; continue;  
	    }
	  }

      	if ( l.in_root && l.in_runs && l.in_run && l.in_dms && depth == 4 &&  xmlStrEqual (name, (const xmlChar *)"dm" ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT: dm.resize(dmSize); l.in_dm = true; l.row_nr = -1; continue; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_dm = false; 
		for(size_t namei=0 ; namei < speciesnames.size() ; namei++ )
		  {   dm.setIdentifier(namei,speciesnames[namei]); }
                return DM_READ;  
	    }
	  }

    	if ( l.in_root && l.in_runs && l.in_run && l.in_identities && depth == 4 && xmlStrEqual (name, (const xmlChar *)"identity" ) &&
               type == XML_READER_TYPE_ELEMENT  ) {
  	        xmlChar * nameStr =  xmlTextReaderGetAttribute(reader,(const xmlChar *)"name" );
                speciesnames.push_back( ( char * ) nameStr );
                xmlFree( nameStr );
	}

      	if ( l.in_root && l.in_runs && l.in_run && depth == 3 &&  xmlStrEqual (name, (const xmlChar *)"dms" ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_dms = true; continue; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_dms = false; continue;  
	    }
	  }

	if (  l.in_root && l.in_runs && depth == 2 && xmlStrEqual (name, (const xmlChar *)"run" ))
	  {
	    if (type == XML_READER_TYPE_ELEMENT ) { l.in_run = true; 
	      xmlChar * dimStr = xmlTextReaderGetAttribute(reader,(const xmlChar *)"dim" );
	      dmSize=atoi((const char *) dimStr );
              xmlFree(dimStr);
	      continue; }; 
            if (type == XML_READER_TYPE_END_ELEMENT ) { l.in_run = false; return END_OF_RUN; }
	  }


      	if ( l.in_root && l.in_runs && l.in_run && depth == 3 &&  xmlStrEqual (name, (const xmlChar *)"identities" ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_identities = true; speciesnames.clear(); continue; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_identities = false; continue;  
	    }
	  }


	if ( l.in_root && depth == 1 &&  xmlStrEqual (name, (const xmlChar *)"runs" ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_runs = true; continue; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_runs = false; return END_OF_RUNS;  
	    }
	  }

	if ( depth == 0 &&  xmlStrEqual (name, (const xmlChar *)"root" ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_root = true; continue; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_root = false; continue; 
	    }
	  }
      }
 }
    return ERROR;
}

