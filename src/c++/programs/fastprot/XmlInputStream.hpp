#ifndef XMLINPUTSTREAM_HPP
#define XMLINPUTSTREAM_HPP

#include <cstdio>

#include <iostream>
#include <fstream>
#include <libxml/xmlreader.h>
#include "DataInputStream.hpp"

#include "../../fileFormatSchema.hpp"

using namespace std;

typedef enum { RUN_NOT_FINISHED = 0, RUN_FINISHED = 1 } streamstatus;

typedef struct {  int in_root; 
  int in_runs;
  int in_run;
  int in_seq; 
 } locator_t;


class XmlInputStream : public DataInputStream
{
public:
   XmlInputStream(char * filename = NULL);
  ~XmlInputStream();

  virtual bool read( std::vector<Sequence> &seqs, std::string & runId, std::vector<std::string> &names, Extrainfos &extrainfos );
  virtual bool readSequences( std::vector<Sequence> &seqs, std::string & runId, Extrainfos &extrainfos );
protected:
  xmlTextReaderPtr reader;
  locator_t l;
  int fd;
};

#endif // XMLINPUTSTREAM_HPP
