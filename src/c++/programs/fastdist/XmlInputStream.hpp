#ifndef XMLINPUTSTREAM_HPP
#define XMLINPUTSTREAM_HPP

#include <cstdio>

#include <iostream>
#include <fstream>
#include <libxml/xmlreader.h>
#include "DataInputStream.hpp"

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
