/*
 * PhylipDmOutputStream.hpp
 *
 *  Created on: Dec 1, 2011
 *      Author: Mehmood Alam Khan
 *      Email: malagori@kth.se
 */

#pragma once

//#include <cstdio>
#include "DataOutputStream.hpp"

class PhylipDmOutputStream : public DataOutputStream
{
public:
  PhylipDmOutputStream(char * filename ) : DataOutputStream(filename) {};
  ~PhylipDmOutputStream() override {};
  void print( StrDblMatrix & dm ) override;
  // changes here for row matrix
  void printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos ) override {};
  void printEndRun() override {};
  void printRow( StrFloRow & dm, std::string name, int row, bool mem_eff_flag) override;
  void printHeader( size_t numNodes ) override;
  void printBootstrapSpliter(size_t numNodes) override;
};

