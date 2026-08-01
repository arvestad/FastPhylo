/*
 * BinaryDmOutputStream.hpp
 *
 *  Created on: Dec 1, 2011
 *      Author: Mehmood Alam Khan
 */

#pragma once

#include "DataOutputStream.hpp"
#include <cstdio>
#include <fstream>


class BinaryDmOutputStream: public DataOutputStream {
public:
	BinaryDmOutputStream(char * filename ) : DataOutputStream(filename) {
	  	if(filename != nullptr) {
	  		writeToCout = false;
	  		fp = nullptr;
	  		delete fp;
	  		file_was_opened = false;
	   	ofs = open_write_binary(filename);
	  	} else {
	  		ofs = &std::cout;
	  		writeToCout = true;
	  	}
		}


	  virtual ~BinaryDmOutputStream() {
	  	if(ofs != nullptr && !writeToCout) {
	  		delete ofs;
	  		ofs = nullptr;
	  	}
		}
	virtual void print( StrDblMatrix & dm ) {};
	  virtual void printHeader( size_t numNodes );
	  virtual void printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos);
	  virtual void printRow( StrFloRow & dm, std::string name, int row, bool mem_eff_flag);
	  virtual void printBootstrapSpliter(size_t numNodes){};
	  virtual void printEndRun() {};

	private:
	  std::ostream *ofs;
	  std::vector<std::string> m_names;
	  bool writeToCout;
};

