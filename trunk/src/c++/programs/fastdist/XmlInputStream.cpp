#include "XmlInputStream.hpp"
#include <cstdio>

#include <libxml/relaxng.h>

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
    filename = "-";
  }
  LIBXML_TEST_VERSION

  reader = xmlReaderForFile(filename,0, XML_PARSE_COMPACT | XML_PARSE_NONET );
  if ( reader == 0 ) { THROW_EXCEPTION("Could not open file"); };

  l.in_root =  false;
  l.in_runs =  false;


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

bool XmlInputStream::read( vector<string> &names, vector<DNA_b128_String> &b128seqs )  
{ 
  std::vector<Sequence> seqs;
  if ( ! readSequences(seqs) ) return false;
  names.clear();names.reserve(seqs.size());
  for( size_t i=0;i<seqs.size();i++) {
    names.push_back(seqs[i].name);
  }
  Sequences2DNA_b128(seqs,b128seqs);
  return true;
}

void 
XmlInputStream::readRunTree( xmlNodePtr tree, std::vector<Sequence> &seqs ) {

  int numSequences = 0;
  unsigned int seqlen;

  xmlNode *node = NULL;
 
  for (node = tree->children; node; node = node->next) {
    if (node->type == XML_ELEMENT_NODE && xmlStrEqual (node->name, (const xmlChar *)"seq" ) ) {
      numSequences++;
    }
  }
  seqs.resize(numSequences);

  int i = 0;
  for (node = tree->children; node; node = node->next) {
    if (node->type == XML_ELEMENT_NODE && xmlStrEqual (node->name, (const xmlChar *)"seq" ) ) {
       Sequence &s = seqs[i];

      xmlChar * name =  xmlGetProp(node, (const xmlChar *) "name") ;
      xmlChar * seq = xmlGetProp(node, (const xmlChar *) "seq") ;

      if ( name == 0 ) THROW_EXCEPTION("failed to read attribute \"name\"");
      if ( seq == 0 ) THROW_EXCEPTION("failed to read attribute \"seq\"");

      s.name =  ( const char *) name;
      s.seq =  ( const char *) seq;

      xmlFree(name);
      xmlFree(seq);
      i++;
    }
  }
}

bool 
XmlInputStream::readSequences( std::vector<Sequence> &seqs ) {
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

	if ( l.in_root && l.in_runs && depth == 2 && xmlStrEqual (name, (const xmlChar *)"run" ))
	  {
	    if ( type == XML_READER_TYPE_ELEMENT ) {

	      xmlNodePtr tree;
	      tree = xmlTextReaderExpand (reader);

              if ( tree == NULL ) { 
                THROW_EXCEPTION("could not expand tree");
                exit(EXIT_FAILURE);
              } 
              readRunTree(tree, seqs ); 

              run_read = true;  continue; 
	    } 
	    if ( type == XML_READER_TYPE_END_ELEMENT ) {
                  return true;
	    } 
	  }

	if ( depth == 0 &&  xmlStrEqual (name, (const xmlChar *)"root" ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_root = true; continue; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_root = false; break; 
	    }
	  }

	if ( l.in_root && depth == 1 &&  xmlStrEqual (name, (const xmlChar *)"runs" ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_runs = true; continue; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_runs = false; continue; 
	    }
	  }
      }
  if (ret == 0 && ! run_read ) {
  THROW_EXCEPTION("YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY");
  exit(EXIT_FAILURE);


    return false;
  }
  THROW_EXCEPTION("failed to parse");
  exit(EXIT_FAILURE);
}
