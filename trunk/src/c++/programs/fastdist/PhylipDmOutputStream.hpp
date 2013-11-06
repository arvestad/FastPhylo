/*
 * PhylipDmOutputStream.hpp
 *
 *  Created on: Dec 1, 2011
 *      Author: Mehmood Alam Khan
 *      Email: malagori@kth.se
 */

#ifndef PHYLIPDMOUTPUTSTREAM_HPP_
#define PHYLIPDMOUTPUTSTREAM_HPP_

//#include <cstdio>
#include "DataOutputStream.hpp"

class PhylipDmOutputStream : public DataOutputStream
{
public:
  PhylipDmOutputStream(char * filename ) : DataOutputStream(filename) {};
  virtual ~PhylipDmOutputStream() {};
  virtual void print( StrDblMatrix & dm );
  // changes here for row matrix
  virtual void printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos ) {};
  virtual void printEndRun() {};
  virtual void printRow( StrFloRow & dm, std::string name, int row, bool mem_eff_flag);
  virtual void printHeader( size_t numNodes );
  virtual void printBootstrapSpliter(size_t numNodes);
};

#endif /* PHYLIPDMOUTPUTSTREAM_HPP_ */
