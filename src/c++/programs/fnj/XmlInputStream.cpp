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

}

void XmlInputStream::read( StrDblMatrix & dm, std::vector<string> & speciesnames, readstatus & status )  
{ 
    const xmlChar *name, *value;

    bool run_read = false;
    int ret;
    status = NOTHING_READ;
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

                for (node = tree->children; node; node = node->next) {
                  if (node->type == XML_ELEMENT_NODE && xmlStrEqual (node->name, (const xmlChar *)"entry" ) ) {
		    speciesnames.push_back(( char * ) node->conent );
                  }
                }
                status = SPECIES_READ; break; // break out of the while loop
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
	    case XML_READER_TYPE_END_ELEMENT:  l.in_runs = false; continue; 
	    }
	  }

	if (  l.in_root && l.in_runs && depth == 2 && xmlStrEqual (name, (const xmlChar *)"run" ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_run = true; continue; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_run = false; continue; 
	    }
	  }
      }

    if (ret == 1 && ( status == SPECIES_READ || status == DM_READ )) {
      return; // successful read
    }
    if (ret == 0 && status == NOTHING_READ  ) { 
      return; // nothing more to read 
    }
 
  THROW_EXCEPTION("failed to parse");
  exit(EXIT_FAILURE);
}


bool 
XmlInputStream::readSpeciesNamesAndOneDM( std::vector<string> & speciesnames, StrDblMatrix & dm ) {

  read( dm, speciesnames, status );

  switch (status){                               
    case SPECIES_READ: { break;  }
    case NOTHING_READ: { return 0; }
    case DM_READ: { 
      THROW_EXCEPTION("failed to parse");
      exit(EXIT_FAILURE);
    }
  }

  std::vector<string> dummy;  

  read( dm, dummy, status );

  switch (status){                               
    case SPECIES_READ: { 
      THROW_EXCEPTION("failed to parse");
      exit(EXIT_FAILURE);
    }
    case NOTHING_READ: { return 0; }
    case DM_READ: { 
      return 1;
    }
  }
}

bool 
XmlInputStream::readSpeciesOneDM( StrDblMatrix & dm ) {

  std::vector<string> dummy;  

  read( dm, dummy, status );

  switch (status){                               
    case SPECIES_READ: { 
      THROW_EXCEPTION("failed to parse");
      exit(EXIT_FAILURE);
    }
    case NOTHING_READ: { return 0; }
    case DM_READ: { 
      return 1;
    }
  }
}
