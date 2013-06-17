/*
 * BinaryDmOutputStream.hpp
 *
 *  Created on: Dec 1, 2011
 *      Author: Mehmood Alam Khan
 */

#ifndef BINARYDMOUTPUTSTREAM_HPP_
#define BINARYDMOUTPUTSTREAM_HPP_

#include "DataOutputStream.hpp"
#include <cstdio>
#include <fstream>


class BinaryDmOutputStream: public DataOutputStream {
public:
	BinaryDmOutputStream(char * filename ) : DataOutputStream(filename) {
	  	if(filename != 0) {
	  		writeToCout = false;
	  		fp = 0;
	  		delete fp;
	  		file_was_opened = false;
	   	ofs = open_write_binary(filename);
	  	} else {
	  		ofs = &std::cout;
	  		writeToCout = true;
	  	}
		}


	  virtual ~BinaryDmOutputStream() {
	  	if(ofs != 0 && !writeToCout) {
	  		delete ofs;
	  		ofs = 0;
	  	}
		}
	virtual void print( StrDblMatrix & dm ) {};
	  virtual void printHeader( size_t numNodes );
	  virtual void printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos);
	  virtual void printRow( StrFloRow & dm, std::string name, int row);
	  virtual void printBootstrapSpliter(size_t numNodes){};
	  virtual void printEndRun() {};

	private:
	  std::ostream *ofs;
	  std::vector<std::string> m_names;
	  bool writeToCout;
};

#endif /* BINARYDMOUTPUTSTREAM_HPP_ */
