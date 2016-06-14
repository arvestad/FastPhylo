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
  readstatus readDM( StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos ) { return ERROR; };
  readstatus readDM( StrFloMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos)  { return ERROR; };
};

/*
  This implementation used to have these two methods as virtual, which is natural with an abstract base class. 
  However, there were warnings for the methods because they had empty bodies ( i.e., their implementation was "{}"). 
  The true way to work with an abstract base class is to have the NULL bodies (i.e., the implemention would be "=0"). 
  But this would lead to compilation errors when subclasses would not implement both methods. 

  My workaround here is to let DataInputStream be a base class and let the default implementations return an error code.

  /arve 2016-06-14
*/

#endif // DATAINPUTSTREAM_HPP
