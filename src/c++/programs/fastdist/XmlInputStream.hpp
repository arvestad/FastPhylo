#ifndef XMLINPUTSTREAM_HPP
#define XMLINPUTSTREAM_HPP

#include <cstdio>

#include <iostream>
#include <fstream>
#include <libxml/xmlreader.h>
#include "DataInputStream.hpp"

static const char relaxngstr_old [] = "<?xml version=\"1.0\"?>\n"
"<element xmlns=\"http://relaxng.org/ns/structure/1.0\" name=\"root\" datatypeLibrary=\"http://www.w3.org/2001/XMLSchema-datatypes\">\n"
" <element name=\"runs\">\n"
"  <zeroOrMore>\n"
"   <element name=\"run\">\n"
"    <attribute name=\"id\">\n"
"     <text/>\n"
"    </attribute>\n"
"    <oneOrMore>\n"
"     <element name=\"seq\">\n"
"      <attribute name=\"seq\">\n"
"       <data type=\"string\">\n"
"        <param name=\"pattern\">[ACGT]+</param>\n"
"       </data>\n"
"      </attribute>\n"
"      <attribute name=\"name\">\n"
"       <text/>\n"
"      </attribute>\n"
"     </element>\n"
"    </oneOrMore>\n"
"   </element>\n"
"  </zeroOrMore>\n"
" </element>\n"
"</element>\n";
  



static const char relaxngstr [] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
"<grammar xmlns=\"http://relaxng.org/ns/structure/1.0\" datatypeLibrary=\"http://www.w3.org/2001/XMLSchema-datatypes\">\n"
"  <start>\n"
"    <element name=\"root\">\n"
"      <element name=\"runs\">\n"
"        <zeroOrMore>\n"
"          <element name=\"run\">\n"
"            <attribute name=\"id\">\n"
"              <text/>\n"
"            </attribute>\n"
"            <oneOrMore>\n"
"              <element name=\"seq\">\n"
"                <attribute name=\"seq\">\n"
"                  <data type=\"string\">\n"
"                    <param name=\"pattern\">[ACGT]+</param>\n"
"                  </data>\n"
"                </attribute>\n"
"                <attribute name=\"name\">\n"
"                  <text/>\n"
"                </attribute>\n"
"                <optional>\n"
"                  <element name=\"extrainfo\">\n"
"                    <ref name=\"anyContent\"/>\n"
"                  </element>\n"
"                </optional>\n"
"              </element>\n"
"            </oneOrMore>\n"
"          </element>\n"
"        </zeroOrMore>\n"
"      </element>\n"
"    </element>\n"
"  </start>\n"
"  <define name=\"anyContent\">\n"
"    <mixed>\n"
"      <zeroOrMore>\n"
"        <choice>\n"
"          <attribute>\n"
"            <anyName/>\n"
"          </attribute>\n"
"          <ref name=\"anyElement\"/>\n"
"        </choice>\n"
"      </zeroOrMore>\n"
"    </mixed>\n"
"  </define>\n"
"  <define name=\"anyElement\">\n"
"    <element>\n"
"      <anyName/>\n"
"      <ref name=\"anyContent\"/>\n"
"    </element>\n"
"  </define>\n"
"</grammar>\n";




typedef enum { RUN_NOT_FINISHED = 0, RUN_FINISHED = 1 } streamstatus;

typedef struct {  int in_root; 
  int in_runs;
  int in_run;
  int in_seq; 
 } locator_t;


class XmlInputStream : public DataInputStream
{
public:
   XmlInputStream(char * filename);
  ~XmlInputStream();
  virtual bool read( std::vector<std::string> &names, std::vector<DNA_b128_String> &b128_strings );
  virtual bool readSequences( std::vector<Sequence> &seqs );
protected:
  xmlTextReaderPtr reader;
  locator_t l;
  int fd;
};

#endif // XMLINPUTSTREAM_HPP
