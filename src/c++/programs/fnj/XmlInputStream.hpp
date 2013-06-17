#ifndef XMLINPUTSTREAM_HPP
#define XMLINPUTSTREAM_HPP

#include <cstdio>
#include <iostream>
#include <fstream>
#include <libxml/xmlreader.h>
#include "DataInputStream.hpp"
#include "fileFormatSchema.hpp"

typedef struct { bool in_root; 
  bool in_runs;
  bool in_run; 
  bool in_identities;
  bool in_identity;
  bool in_dms; 
  bool in_dm;
  bool in_row; 
  int row_nr;
  int entry_nr; 
 } locator_t;

class XmlInputStream : public DataInputStream {
public:
  XmlInputStream(char *filename);
  ~XmlInputStream();
  readstatus readDM( StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos);
protected:
  xmlTextReaderPtr reader;
  locator_t l;
  int fd;
  int dmSize;
};

#endif // XMLINPUTSTREAM_HPP
