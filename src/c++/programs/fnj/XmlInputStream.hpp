#ifndef XMLINPUTSTREAM_HPP
#define XMLINPUTSTREAM_HPP

#include <cstdio>

#include <iostream>
#include <fstream>
#include <libxml/xmlreader.h>
#include "DataInputStream.hpp"

static const char relaxngstr [] = "<?xml version=\"1.0\"?>\n"
"<element xmlns=\"http://relaxng.org/ns/structure/1.0\" name=\"root\" datatypeLibrary=\"http://www.w3.org/2001/XMLSchema-datatypes\">\n"
" <element name=\"runs\">\n"
"  <zeroOrMore>\n"
"   <element name=\"run\">\n"
"    <attribute name=\"dim\">\n"
"     <data type=\"integer\"/>\n"
"    </attribute>\n"
"    <element name=\"identities\">\n"
"     <oneOrMore>\n"
"      <element name=\"identity\">\n"
"       <attribute name=\"name\">\n"
"        <text/>\n"
"       </attribute>\n"
"      </element>\n"
"     </oneOrMore>\n"
"    </element>\n"
"    <element name=\"dms\">\n"
"     <oneOrMore>\n"
"      <element name=\"dm\">\n"
"       <oneOrMore>\n"
"        <element name=\"row\">\n"
"         <oneOrMore>\n"
"          <element name=\"entry\">\n"
"           <data type=\"float\"/>\n"
"          </element>\n"
"         </oneOrMore>\n"
"        </element>\n"
"       </oneOrMore>\n"
"      </element>\n"
"     </oneOrMore>\n"
"    </element>\n"
"   </element>\n"
"  </zeroOrMore>\n"
" </element>\n"
"</element>\n";


typedef struct { bool in_root; 
  bool in_runs;
  bool in_run; 
  bool in_identities;
  bool in_dms; 
  bool in_dm;
  bool in_row; 
  int row_nr;
  int entry_nr; 
 } locator_t;


class XmlInputStream : public DataInputStream
{
public:
   XmlInputStream(char * filename);
  ~XmlInputStream();

virtual readstatus readDM( StrDblMatrix & dm );

protected:

  xmlTextReaderPtr reader;
  locator_t l;
  int fd;
  int dmSize;

  std::vector<std::string> speciesnames;
};

#endif // XMLINPUTSTREAM_HPP
