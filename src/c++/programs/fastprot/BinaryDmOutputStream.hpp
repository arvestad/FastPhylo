/*
 * BinaryDmOutputStream.hpp
 *
 *  Created on: March 18, 2013
 *      Author: Henric Zazzi
 */

#pragma once

#include "DataOutputStream.hpp"
#include <cstdio>
#include <fstream>
#include <memory>

using namespace std;

class BinaryDmOutputStream: public DataOutputStream {
public:
	BinaryDmOutputStream(char * filename);
	void printHeader( size_t numNodes );
	void printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos);
	void print(StrDblMatrix &dm);

private:
	// Owns the file when writing to disk; unused (writeToCout) when
	// writing to stdout, in which case ofs aliases &cout instead.
	std::unique_ptr<std::ofstream> file_stream;
	ostream *ofs;
	vector<string> m_names;
	bool writeToCout;
};

