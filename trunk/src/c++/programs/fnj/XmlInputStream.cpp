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

bool XmlInputStream::read( StrDblMatrix & dm,  str2int_hashmap & name2id )  
{ 
    const xmlChar *name, *value;

    bool run_read = false;
    int ret;
    while ( ( ret = xmlTextReaderRead(reader)) == 1 )
      {
	int depth = xmlTextReaderDepth(reader);
	int type = xmlTextReaderNodeType(reader);
	name = xmlTextReaderConstName(reader);

	if ( l.in_root && l.in_runs && depth == 2 && xmlStrEqual (name, (const xmlChar *)"run" ))
	  {
	    if ( type == XML_READER_TYPE_ELEMENT ) {

	      xmlNodePtr tree;
	      tree = xmlTextReaderExpand (reader);
              readRun(tree, dm, name2id  ); 
 
              run_read = true;  break; // break out of the while loop
	    } 
	  }

	if ( depth == 0 &&  xmlStrEqual (name, (const xmlChar *)"root" ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_root = true; break; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_root = false; break; 
	    }
	  }

	if ( l.in_root && depth == 1 &&  xmlStrEqual (name, (const xmlChar *)"runs" ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_runs = true; break; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_runs = false; break; 
	    }
	  }
      }

  if (ret == 1 && run_read ) {
    return true;
  }

  if (ret == 0 && ! run_read ) {
    return false;
  }
  THROW_EXCEPTION("failed to parse");
  exit(EXIT_FAILURE);
}

void 
XmlInputStream::readRun( xmlNodePtr tree,  StrDblMatrix & dm,  str2int_hashmap & name2id  ) {

  int numSequences = 0;
  unsigned int seqlen;

  xmlNode *node = NULL;
 
  for (node = tree->children; node; node = node->next) {
    if (node->type == XML_ELEMENT_NODE && xmlStrEqual (node->name, (const xmlChar *)"seq" ) ) {
      numSequences++;
    }
  }
  //  seqs.resize(numSequences);

  int i = 0;
  for (node = tree->children; node; node = node->next) {
    if (node->type == XML_ELEMENT_NODE && xmlStrEqual (node->name, (const xmlChar *)"seq" ) ) {
      //       Sequence &s = seqs[i];

      xmlChar * name =  xmlGetProp(node, (const xmlChar *) "name") ;
      xmlChar * seq = xmlGetProp(node, (const xmlChar *) "seq") ;

      if ( name == 0 ) THROW_EXCEPTION("failed to read attribute \"name\"");
      if ( seq == 0 ) THROW_EXCEPTION("failed to read attribute \"seq\"");

      //      s.name =  ( const char *) name;
      //      s.seq =  ( const char *) seq;

      xmlFree(name);
      xmlFree(seq);
      i++;
    }
  }
}
