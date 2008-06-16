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

  dmSize = 0;
}

readstatus  
XmlInputStream::readDM( StrDblMatrix & dm ) {
{ 
    const xmlChar *name, *value;

    bool run_read = false;
    int ret;

    while ( ( ret = xmlTextReaderRead(reader)) == 1 )
      {
	int depth = xmlTextReaderDepth(reader);
	int type = xmlTextReaderNodeType(reader);
	name = xmlTextReaderConstName(reader);

	if ( l.in_root && l.in_runs && l.in_run && depth == 3 )
	  {
	    if ( xmlStrEqual (name, (const xmlChar *)"species" ) ) {
	      if ( type == XML_READER_TYPE_ELEMENT ) {
		xmlNodePtr tree;
		tree = xmlTextReaderExpand (reader);

                xmlNode *node = NULL;
		int i=0;
                speciesnames.clear();
                for (node = tree->children; node; node = node->next) {
                  if (node->type == XML_ELEMENT_NODE && xmlStrEqual (node->name, (const xmlChar *)"entry" ) ) {
		    char * species =  ( char * ) xmlNodeGetContent( node );
    		    speciesnames.push_back( species );
                    i++;
                  }
                }
	      } 
	    } 
	    if ( xmlStrEqual (name, (const xmlChar *)"dm" ) ) {
	      if ( type == XML_READER_TYPE_ELEMENT ) {
		xmlNodePtr tree;
		tree = xmlTextReaderExpand (reader);
                dm.resize(dmSize);

                xmlNode *node = NULL;
                int rowIter = 0;
                for (node = tree->children; node; node = node->next) {
                  if (node->type == XML_ELEMENT_NODE && xmlStrEqual (node->name, (const xmlChar *)"row" ) ) {
		    //		    speciesnames.push_back(( char * ) node->content );
                    int entryIter = 0; 
                    for (xmlNode *node2  = node->children; node2; node2 = node2->next) {
                      if (node2->type == XML_ELEMENT_NODE && xmlStrEqual (node2->name, (const xmlChar *)"entry" ) ) {
                	//	dm.setDistance(rowIter,entryIter,  atof(( char * ) node2->content ));
			entryIter++;
		      }
		    }
                    rowIter++;
                  }
                }
		for(size_t namei=0 ; namei < speciesnames.size() ; namei++ )
		  {	    dm.setIdentifier(namei,speciesnames[namei]); }

                return DM_READ; // break out of the while loop
	      } 
	    } 
	  }

	if ( depth == 0 &&  xmlStrEqual (name, (const xmlChar *)"root" ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_root = true; continue; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_root = false; continue; 
	    }
	  }

	if ( l.in_root && depth == 1 &&  xmlStrEqual (name, (const xmlChar *)"runs" ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_runs = true; continue; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_runs = false; return END_OF_RUNS;  
	    }
	  }

	if (  l.in_root && l.in_runs && depth == 2 && xmlStrEqual (name, (const xmlChar *)"run" ))
	  {
	    if (type == XML_READER_TYPE_ELEMENT ) { l.in_run = true; 
	      char * dimStr = ( char * ) xmlTextReaderGetAttribute(reader,(const xmlChar *)"dim" );
	      dmSize=atoi((const char *) dimStr );
	      continue; }; 
            if (type == XML_READER_TYPE_END_ELEMENT ) { l.in_run = false; return END_OF_RUN; }
	  }
      }
 }
    return ERROR;
}

