#pragma once

#include <cstdio>
#include "fastphylo/dna/Sequences2DistanceMatrix.hpp"
#include "fastphylo/core/Exception.hpp"
#include "fastphylo/io/Extrainfos.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

enum readstatus { DM_READ = 1, END_OF_RUN = 2, END_OF_RUNS = 3, ERROR = 4 };

class DataInputStream
{
public:
  virtual ~DataInputStream() {};
  virtual readstatus readDM(StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos ) = 0;
};

/*
  This implementation used to have readDM as virtual with an empty method
  body (i.e., its implementation was "{}"). That created warnings/errors
  with more modern compilers. The true way to work with an abstract base
  class is to have a NULL body (i.e., "=0"), which is the current solution.

  /arve 2016-06-14
*/

