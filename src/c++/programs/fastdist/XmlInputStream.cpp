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

XmlInputStream::XmlInputStream(char * filename)  
{
  //
  //  Re: [xml] Why does "-" read from stdin?
  //  http://mail.gnome.org/archives/xml/2007-February/msg00005.html
  //
  // Comment: The special treatment of the filename "-" in the libxml api, is not  
  // a good designed api, but now when it is there let us use it.

  if (filename && !strncmp(filename,"-",1)) {
    THROW_EXCEPTION("file name \"-\" is not allowed. See \n http://mail.gnome.org/archives/xml/2007-February/msg00005.html \n  ");
  }

  if ( filename == nullptr ) {
    filename = const_cast<char *>("-");
  }
  LIBXML_TEST_VERSION

  reader = xmlReaderForFile(filename,nullptr, XML_PARSE_COMPACT | XML_PARSE_NONET );
  if ( reader == nullptr ) { THROW_EXCEPTION("Could not open file"); };

  l.in_root = false;
  l.in_runs = false;
  l.in_run  = false;
  l.in_seq  = false;

  xmlRelaxNGParserCtxtPtr parserctxt;
  size_t len = strlen(fastphylo_sequence_xml_relaxngstr);
  parserctxt = xmlRelaxNGNewMemParserCtxt(fastphylo_sequence_xml_relaxngstr,len);
  xmlRelaxNGSetParserErrors(parserctxt,reinterpret_cast<xmlRelaxNGValidityErrorFunc>(fprintf), reinterpret_cast<xmlRelaxNGValidityWarningFunc>(fprintf), stderr);
  xmlRelaxNGPtr schema = nullptr;
  schema = xmlRelaxNGParse(parserctxt);
  xmlRelaxNGFreeParserCtxt(parserctxt);

  if ( xmlTextReaderRelaxNGSetSchema( reader,  schema ) != 0 )  { 
    THROW_EXCEPTION("failed to set relax ng schema");
    exit(EXIT_FAILURE);
  }
}

bool XmlInputStream::read(  std::vector<DNA_b128_String> &b128seqs, std::string & runId, std::vector<std::string> &names, Extrainfos &extrainfos )  
{ 
  std::vector<Sequence> seqs;
  if ( ! readSequences(seqs, runId, extrainfos) ) return false;
  names.clear();names.reserve(seqs.size());
  for( size_t i=0;i<seqs.size();i++) {
    names.push_back(seqs[i].name);
  }
  Sequences2DNA_b128(seqs,b128seqs);
  return true;
}

bool
XmlInputStream::readSequences( std::vector<Sequence> &seqs, std::string & runId, Extrainfos &extrainfos  ) {
    const xmlChar *name, *value;

    bool run_read = false;
    int ret;
    int numSequences = 0;
    while ( ( ret = xmlTextReaderRead(reader)) == 1 )
      {
 	if ( xmlTextReaderIsValid( reader ) != 1 ) { 
          THROW_EXCEPTION("xml input does not validate");
          exit(EXIT_FAILURE);
        } 

	int depth = xmlTextReaderDepth(reader);
	int type = xmlTextReaderNodeType(reader);
	name = xmlTextReaderConstName(reader);

	if ( l.in_root && l.in_runs && l.in_run && depth == 3 && xmlStrEqual (name, reinterpret_cast<const xmlChar *>("seq") ))
	  {
	    if ( type == XML_READER_TYPE_ELEMENT ) {
              l.in_seq = true;
	      numSequences++;
              seqs.resize(numSequences);
	      extrainfos.push_back( std::string() );
	      //              extrainfos.resize(numSequences);
	      //    extrainfos[numSequences-1]=NULL;
              Sequence &s = seqs[numSequences-1];
 
              xmlChar * name =  xmlTextReaderGetAttribute(reader, reinterpret_cast<const xmlChar *>("name")) ;
              xmlChar * seq = xmlTextReaderGetAttribute(reader, reinterpret_cast<const xmlChar *>("seq")) ;

              if ( name == nullptr ) THROW_EXCEPTION("failed to read attribute \"name\"");
              if ( seq == nullptr ) THROW_EXCEPTION("failed to read attribute \"seq\"");

              s.name =  reinterpret_cast<const char *>(name);
              s.seq = reinterpret_cast<const char *>(seq);
              xmlFree(name);
              xmlFree(seq);
              continue;
	    } 
	    if ( type == XML_READER_TYPE_END_ELEMENT ) {
	          l.in_seq = false;
                  continue;
	    } 
        }

	if ( l.in_root && l.in_runs && l.in_run && l.in_seq && depth == 4 && xmlStrEqual (name, reinterpret_cast<const xmlChar *>("extrainfo") ))
	  {
	    if ( type == XML_READER_TYPE_ELEMENT ) {
 
	      //              xmlNodePtr extrainfoPtr = xmlTextReaderExpand(reader);
	      // if ( extrainfoPtr == 0 ) THROW_EXCEPTION("could not xmlTextReaderExpand \"extrainfo\"");
              xmlChar * outerStr = xmlTextReaderReadOuterXml(reader);

	      //              extrainfos[numSequences-1]=outerStr;
              extrainfos.back()= reinterpret_cast<char *>(outerStr);
              xmlFree(outerStr);

	      //              xmlChar * extrainfo = extrainfos.back();
	      //            extrainfo = ;
	      //  char * m = ( char *) outerStr;
	      // std::cout <<  m   << std::endl;

			    //              xmlFree(outerStr);
              continue;

	    } 
	    if ( type == XML_READER_TYPE_END_ELEMENT ) {
                  continue;
	    } 
        }

	if ( depth == 0 &&  xmlStrEqual (name, reinterpret_cast<const xmlChar *>("root") ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_root = true; continue; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_root = false; break; 
	    }
	  }

	if ( l.in_root && depth == 1 &&  xmlStrEqual (name, reinterpret_cast<const xmlChar *>("runs") ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  l.in_runs = true; continue; 
	    case XML_READER_TYPE_END_ELEMENT:  l.in_runs = false; continue; 
	    }
	  }

	if ( l.in_root && l.in_runs && depth == 2 && xmlStrEqual (name, reinterpret_cast<const xmlChar *>("run") ))
	  {
	    switch (type) {
	    case XML_READER_TYPE_ELEMENT:  {
	      l.in_run = true;
	      extrainfos.clear(); 
              xmlChar * id =  xmlTextReaderGetAttribute(reader, reinterpret_cast<const xmlChar *>("id")) ;
              if ( id == nullptr ) THROW_EXCEPTION("failed to read attribute \"id\"");
              runId = reinterpret_cast<const char *>(id);
              xmlFree(id); continue;
	    }
	    case XML_READER_TYPE_END_ELEMENT:  {
	      l.in_run = false; 
	      return true;
	    }
            }
	  }
      }
    if (ret == 0 ) {
      if ( l.in_root ) {
        THROW_EXCEPTION("failed to parse");
        exit(EXIT_FAILURE);
      } 
      else { 
        return false;
      }
    }
    return false;		// It cannot be good that we reach here. /arve
}
