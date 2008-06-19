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
"    <attribute name=\"id\">\n"
"     <text/>\n"
"    </attribute>\n"
"    <oneOrMore>\n"
"     <element name=\"seq\">\n"
"      <attribute name=\"seq\">\n"
"       <data type=\"string\">"
"        <param name=\"pattern\">[ACGT]+</param>"
"       </data>"
"      </attribute>\n"
"      <attribute name=\"name\">\n"
"       <text/>\n"
"      </attribute>\n"
"     </element>\n"
"    </oneOrMore>\n"
"   </element>\n"
"  </zeroOrMore>\n"
" </element>\n"
" </element>\n";
  



typedef enum { RUN_NOT_FINISHED = 0, RUN_FINISHED = 1 } streamstatus;

typedef struct {  int in_root ; 
  int in_runs ; 
 } locator_t;


class XmlInputStream : public DataInputStream
{
public:
   XmlInputStream(char * filename);
  ~XmlInputStream();
  bool read( std::vector<std::string> &names, std::vector<DNA_b128_String> &b128_strings );
  bool readSequences( std::vector<Sequence> &seqs );
protected:
  void readRunTree( xmlNodePtr tree, std::vector<Sequence> &seqs );
  xmlTextReaderPtr reader;
  locator_t l;
  int fd;
};

#endif // XMLINPUTSTREAM_HPP
