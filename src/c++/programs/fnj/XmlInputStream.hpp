#ifndef XMLINPUTSTREAM_HPP
#define XMLINPUTSTREAM_HPP

#include <cstdio>

#include <iostream>
#include <fstream>
#include <libxml/xmlreader.h>
#include "DataInputStream.hpp"



typedef struct {  int in_root ; 
  int in_runs ;
  int in_run ; 
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
