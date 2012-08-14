#ifndef DATAINPUTSTREAM_HPP
#define DATAINPUTSTREAM_HPP

#include <cstdio>
#include "Sequences2DistanceMatrix.hpp"
#include "Exception.hpp"
#include "Extrainfos.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <libxml/tree.h>

typedef enum { DM_READ = 1, END_OF_RUN = 2, END_OF_RUNS = 3, ERROR = 4 } readstatus;

class DataInputStream
{
public:
  DataInputStream();
  virtual ~DataInputStream() {};
  virtual readstatus readDM( StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos )=0 ;
  virtual readstatus readFloatDM( StrFloMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos )=0 ;

};

/*class PhylipDmInputStream : public DataInputStream
{
public:
  PhylipDmInputStream(char * filename );
  ~PhylipDmInputStream();
  virtual readstatus readDM( StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos );

protected:
  std::istream * fp;
  std::ifstream fin;
  bool file_was_opened;

 
};
*/

#endif // DATAINPUTSTREAM_HPP
