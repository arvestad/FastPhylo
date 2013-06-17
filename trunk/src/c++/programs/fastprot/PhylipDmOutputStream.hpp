#ifndef PHYLIPDMOUTPUTSTREAM_HPP
#define PHYLIPDMOUTPUTSTREAM_HPP

#include <cstdio>
#include <vector>
#include <string>
#include <math.h>
#include "Exception.hpp"
#include "DistanceMatrix.hpp"
#include "DataOutputStream.hpp"
#include "DistanceRow.hpp"
#include "Extrainfos.hpp"

class PhylipDmOutputStream : public DataOutputStream {
public:
  PhylipDmOutputStream(char * filename ) : DataOutputStream(filename) {};
  void print( StrDblMatrix & dm );
  void printSD( StrDblMatrix & dm );
  static void printPHYLIPfastSD(const StrDblMatrix &dm, FILE *out, bool writeXml, bool writeXmlSD );
};

#endif // PHYLIPDMOUTPUTSTREAM_HPP
