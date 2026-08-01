#pragma once

#include <cstdio>
#include <vector>
#include <string>
#include <math.h>
#include "fastphylo/core/Exception.hpp"
#include "fastphylo/core/DistanceMatrix.hpp"
#include "DataOutputStream.hpp"
#include "fastphylo/core/DistanceRow.hpp"
#include "Extrainfos.hpp"

class PhylipDmOutputStream : public DataOutputStream {
public:
  PhylipDmOutputStream(char * filename ) : DataOutputStream(filename) {};
  void print( StrDblMatrix & dm ) override;
  void printSD( StrDblMatrix & dm ) override;
  static void printPHYLIPfastSD(const StrDblMatrix &dm, FILE *out, bool writeXml, bool writeXmlSD );
};

