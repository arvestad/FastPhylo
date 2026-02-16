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
  virtual ~DataInputStream() {};
  virtual readstatus readDM(StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos ) = 0;
//  virtual readstatus readDM(StrFloMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos) = 0;
};

/*
  This implementation used to have these two methods as virtual, but with empty method bodies ( i.e., their implementation was "{}").
  That created warnings/errors with more modern compilers. The true way to work with an abstract base class is to have the NULL bodies (i.e., the implemention would be "=0"), which is the current solution. However, to make that work, I have added empty implementations 
  (although with error messages, should they be called) for the one method per subclass. This is because it is only in
  BinaryInputStream that we want StrFloMatrix for the first argument.

  /arve 2016-06-14
*/

#endif // DATAINPUTSTREAM_HPP
