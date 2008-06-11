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
  bool read( StrDblMatrix & dm,  str2int_hashmap & name2id );
protected:
  void readRun( xmlNodePtr tree,  StrDblMatrix & dm,  str2int_hashmap & name2id );
  xmlTextReaderPtr reader;
  locator_t l;
  int fd;
};

#endif // XMLINPUTSTREAM_HPP
