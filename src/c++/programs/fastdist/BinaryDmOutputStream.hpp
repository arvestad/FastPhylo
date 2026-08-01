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
#include <memory>


class BinaryDmOutputStream: public DataOutputStream {
public:
	BinaryDmOutputStream(char * filename ) : DataOutputStream(filename) {
	  	if(filename != nullptr) {
	  		writeToCout = false;
	  		// The base DataOutputStream(filename) constructor already opened fp
	  		// via fopen() (open_write_file()); we don't need it (we write
	  		// through ofs instead), so close it properly. This used to be
	  		// `fp = nullptr; delete fp;` - nulling the pointer and then
	  		// deleting the null (a harmless no-op) without ever closing the
	  		// file fopen()'d by the base constructor, silently leaking that
	  		// file descriptor on every invocation of `-O binary -o <file>`.
	  		fclose(fp);
	  		fp = nullptr;
	  		file_was_opened = false;
	  		file_stream.reset(open_write_binary(filename));
	  		ofs = file_stream.get();
	  	} else {
	  		ofs = &std::cout;
	  		writeToCout = true;
	  	}
		}

	virtual void print( StrDblMatrix & dm ) {};
	  virtual void printHeader( size_t numNodes );
	  virtual void printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos);
	  virtual void printRow( StrFloRow & dm, std::string name, int row, bool mem_eff_flag);
	  virtual void printBootstrapSpliter(size_t numNodes){};
	  virtual void printEndRun() {};

	private:
	  // Owns the file when writing to disk; unused (writeToCout) when
	  // writing to stdout, in which case ofs aliases &cout instead.
	  std::unique_ptr<std::ofstream> file_stream;
	  std::ostream *ofs;
	  std::vector<std::string> m_names;
	  bool writeToCout;
};

