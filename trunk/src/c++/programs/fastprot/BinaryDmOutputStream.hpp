/*
 * BinaryDmOutputStream.hpp
 *
 *  Created on: March 18, 2013
 *      Author: Henric Zazzi
 */

#ifndef BINARYDMOUTPUTSTREAM_HPP_
#define BINARYDMOUTPUTSTREAM_HPP_

#include "DataOutputStream.hpp"
#include <cstdio>
#include <fstream>

using namespace std;

class BinaryDmOutputStream: public DataOutputStream {
public:
	BinaryDmOutputStream(char * filename);
	~BinaryDmOutputStream();
	void printHeader( size_t numNodes );
	void printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos);
	void print(StrDblMatrix &dm);

private:
	ostream *ofs;
	vector<string> m_names;
	bool writeToCout;
};

#endif /* BINARYDMOUTPUTSTREAM_HPP_ */
