#include "XmlInputStream.hpp"
#include <cstdio>

using namespace std;

XmlInputStream::~XmlInputStream() {
  if (reader)
    xmlFreeTextReader(reader);
  xmlCleanupParser();
}

XmlInputStream::XmlInputStream(char * filename) {
  //
  //  Re: [xml] Why does "-" read from stdin?
  //  http://mail.gnome.org/archives/xml/2007-February/msg00005.html
  //
  // Comment: The special treatment of the filename "-" in the libxml api, is not  
  // a good designed api, but now when it is there let us use it.

  if ( !strncmp(filename, "-",1) ) {
    THROW_EXCEPTION("file name \"-\" is not allowed. See \n http://mail.gnome.org/archives/xml/2007-February/msg00005.html \n  ");
  }
  if ( filename==NULL) {
    filename = (char *) "-";
  }
  LIBXML_TEST_VERSION
  reader = xmlReaderForFile(filename,0, XML_PARSE_COMPACT | XML_PARSE_NONET );
  if ( reader == 0 )
    THROW_EXCEPTION("Could not open file");
  l.in_root =  false;
  l.in_runs =  false;
  l.in_run = false;
  l.in_identities = false;
  l.in_identity = false;
  l.row_nr = -1;
  l.entry_nr = -1;
  dmSize = 0;
  xmlRelaxNGParserCtxtPtr parserctxt;
  size_t len = strlen(fastphylo_distance_matrix_xml_relaxngstr);
  parserctxt = xmlRelaxNGNewMemParserCtxt(fastphylo_distance_matrix_xml_relaxngstr,len);
  xmlRelaxNGSetParserErrors(parserctxt,(xmlRelaxNGValidityErrorFunc) fprintf, (xmlRelaxNGValidityWarningFunc) fprintf, stderr);
  xmlRelaxNGPtr schema = NULL;
  schema = xmlRelaxNGParse(parserctxt);
  xmlRelaxNGFreeParserCtxt(parserctxt);
  if ( xmlTextReaderRelaxNGSetSchema( reader,  schema ) != 0 )  { 
    THROW_EXCEPTION("failed to set relax ng schema");
    exit(EXIT_FAILURE);
  }
}

readstatus  XmlInputStream::readDM( StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos ) {
  const xmlChar *name, *value;
  bool run_read = false;
  int ret;
  int nr_of_ids = 0;

  while ( ( ret = xmlTextReaderRead(reader)) == 1 ) {
    if ( xmlTextReaderIsValid( reader ) != 1 ) {
      THROW_EXCEPTION("xml input does not validate");
      exit(EXIT_FAILURE);
    }
    int depth = xmlTextReaderDepth(reader);
    int type = xmlTextReaderNodeType(reader);
    name = xmlTextReaderConstName(reader);
    if (l.in_root && l.in_runs && l.in_run && l.in_dms && l.in_dm &&  l.in_row &&
        depth == 6 &&  xmlStrEqual (name, (const xmlChar *)"entry" )  && type == XML_READER_TYPE_ELEMENT) {
      l.entry_nr++;
      xmlChar * distanceStr = xmlTextReaderReadString(reader);
      float distance =  atof( ( char * ) distanceStr );
      dm.setDistance(l.row_nr,l.entry_nr, distance );
      xmlFree(distanceStr);
      continue;
    }
    if (l.in_root && l.in_runs && l.in_run && l.in_dms && l.in_dm &&
         depth == 5 &&  xmlStrEqual (name, (const xmlChar *)"row" )) {
	    switch (type) {
				case XML_READER_TYPE_ELEMENT:
					l.in_row = true;
					l.row_nr++;
					l.entry_nr = -1;
					continue;
				case XML_READER_TYPE_END_ELEMENT:
					l.in_row = false;
					continue;
	    }
	  }
		if ( l.in_root && l.in_runs && l.in_run && l.in_dms && depth == 4 &&
				 xmlStrEqual (name, (const xmlChar *)"dm" )) {
	    switch (type) {
				case XML_READER_TYPE_ELEMENT:
					dm.resize(dmSize);
					l.in_dm = true;
					l.row_nr = -1;
					continue;
				case XML_READER_TYPE_END_ELEMENT:
					l.in_dm = false;
					for(size_t namei=0; namei < names.size(); namei++)
						dm.setIdentifier(namei,names[namei]);
					return DM_READ;
	    }
	  }
		if ( l.in_root && l.in_runs && l.in_run && l.in_identities &&
				 depth == 4 && xmlStrEqual (name, (const xmlChar *)"identity" ) ) {
			switch (type) {
				case XML_READER_TYPE_ELEMENT:
					nr_of_ids++;
					l.in_identity = true;
					extrainfos.push_back( std::string() );
					//                extrainfos.resize(nr_of_ids);
					//                extrainfos[nr_of_ids-1]=NULL;
					xmlChar *nameStr;
					nameStr=xmlTextReaderGetAttribute(reader,(const xmlChar *)"name" );
					names.push_back( (char *) nameStr );
					xmlFree(nameStr);
					continue;
				case XML_READER_TYPE_END_ELEMENT:
					l.in_identity = false;
					continue;
			}
		}
		if ( l.in_root && l.in_runs && l.in_run && l.in_identities && l.in_identity &&
				 depth == 5 && xmlStrEqual (name, (const xmlChar *)"extrainfo" )) {
	    switch (type) {
				case XML_READER_TYPE_ELEMENT:
					xmlChar *outerStr;
					outerStr = xmlTextReaderReadOuterXml(reader);
					extrainfos.back() = ( char * ) outerStr;
					xmlFree(outerStr);
					continue;
				case XML_READER_TYPE_END_ELEMENT:
					continue;
	    }
		}
		if ( l.in_root && l.in_runs && l.in_run && depth == 3 &&
				 xmlStrEqual (name, (const xmlChar *)"dms" )) {
	    switch (type) {
				case XML_READER_TYPE_ELEMENT:
					l.in_dms = true;
					continue;
				case XML_READER_TYPE_END_ELEMENT:
					l.in_dms = false;
					continue;
	    }
	  }
		if (  l.in_root && l.in_runs && depth == 2 && xmlStrEqual (name, (const xmlChar *)"run" )) {
			if (type == XML_READER_TYPE_ELEMENT ) {
				l.in_run = true;
	      xmlChar * dimStr = xmlTextReaderGetAttribute(reader,(const xmlChar *)"dim" );
	      xmlChar * idStr = xmlTextReaderGetAttribute(reader,(const xmlChar *)"id" );
				if ( dimStr == 0 )
					THROW_EXCEPTION("failed to read attribute \"dim\"");
				if ( idStr == 0 )
					THROW_EXCEPTION("failed to read attribute \"id\"");
				dmSize=atoi((const char *) dimStr);
				runId = (const char *) idStr;
				xmlFree(idStr);
				xmlFree(dimStr);
				continue;
			}
			if (type == XML_READER_TYPE_END_ELEMENT ) {
				l.in_run = false;
				return END_OF_RUN;
			}
	  }
		if ( l.in_root && l.in_runs && l.in_run && depth == 3 &&
				 xmlStrEqual (name, (const xmlChar *)"identities" )) {
	    switch (type) {
				case XML_READER_TYPE_ELEMENT:
					l.in_identities = true;
					names.clear();
					extrainfos.clear();
					continue;
				case XML_READER_TYPE_END_ELEMENT:
					l.in_identities = false;
					continue;
	    }
	  }
		if ( l.in_root && depth == 1 &&  xmlStrEqual (name, (const xmlChar *)"runs" )) {
	    switch (type) {
				case XML_READER_TYPE_ELEMENT:
					l.in_runs = true;
					continue;
				case XML_READER_TYPE_END_ELEMENT:
					l.in_runs = false;
					return END_OF_RUNS;
	    }
	  }
		if ( depth == 0 &&  xmlStrEqual (name, (const xmlChar *)"root" )) {
	    switch (type) {
				case XML_READER_TYPE_ELEMENT:
					l.in_root = true;
					continue;
				case XML_READER_TYPE_END_ELEMENT:
					l.in_root = false;
					continue;
	    }
	  }
	}
	return ERROR;
}


readstatus  XmlInputStream::readDM( StrFloMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos ) {
  std::cerr << "XmlInputStream::readDM(StrFloMatrix, ...) -- Not implemented!" << endl;
  std::exit(-1);
}
