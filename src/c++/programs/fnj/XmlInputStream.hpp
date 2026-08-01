#pragma once

#include <cstdio>
#include <iostream>
#include <fstream>
#include <libxml/xmlreader.h>
#include "DataInputStream.hpp"
#include "fileFormatSchema.hpp"

struct locator_t {
  bool in_root;
  bool in_runs;
  bool in_run;
  bool in_identities;
  bool in_identity;
  bool in_dms;
  bool in_dm;
  bool in_row;
  int row_nr;
  int entry_nr;
};

class XmlInputStream : public DataInputStream {
public:
  XmlInputStream(char *filename);
  ~XmlInputStream() override;
  readstatus readDM( StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos) override;
  // Not an override - DataInputStream only declares the StrDblMatrix
  // overload as the (pure) virtual interface; this StrFloMatrix
  // overload is specific to this class (see DataInputStream.hpp's
  // 2016-06-14 comment on why the two aren't symmetric).
  readstatus readDM( StrFloMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos);
protected:
  xmlTextReaderPtr reader;
  locator_t l;
  int fd;
  int dmSize;
};

