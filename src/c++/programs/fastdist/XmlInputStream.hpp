#pragma once

#include <cstdio>

#include <iostream>
#include <fstream>
#include <libxml/xmlreader.h>
#include "DataInputStream.hpp"

#include "fileFormatSchema.hpp"




enum streamstatus { RUN_NOT_FINISHED = 0, RUN_FINISHED = 1 };

struct locator_t {
  int in_root;
  int in_runs;
  int in_run;
  int in_seq;
};


class XmlInputStream : public DataInputStream
{
public:
   XmlInputStream(char * filename = nullptr);
  ~XmlInputStream();

  virtual bool read( std::vector<DNA_b128_String> &b128_strings, std::string & runId, std::vector<std::string> &names, Extrainfos &extrainfos );
  virtual bool readSequences( std::vector<Sequence> &seqs, std::string & runId, Extrainfos &extrainfos );
protected:
  xmlTextReaderPtr reader;
  locator_t l;
  int fd;
};

