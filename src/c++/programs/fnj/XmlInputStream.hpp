#ifndef XMLINPUTSTREAM_HPP
#define XMLINPUTSTREAM_HPP

#include <cstdio>

#include <iostream>
#include <fstream>
#include <libxml/xmlreader.h>
#include "DataInputStream.hpp"

typedef enum { SPECIES_READ = 0, DM_READ = 1, NOTHING_READ = 3 } readstatus;

typedef struct {  int in_root ; 
  int in_runs ;
  int in_run ; 
 } locator_t;


class XmlInputStream : public DataInputStream
{
public:
   XmlInputStream(char * filename);
  ~XmlInputStream();

  bool readSpeciesNamesAndOneDM( std::vector<std::string> & speciesnames, StrDblMatrix & dm );
  bool readSpeciesOneDM( StrDblMatrix & dm );

protected:
  void read( StrDblMatrix & dm, std::vector<std::string> & speciesnames, readstatus & status );

  xmlTextReaderPtr reader;
  locator_t l;
  int fd;
};

#endif // XMLINPUTSTREAM_HPP
