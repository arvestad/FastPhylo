/*
 * BinaryDmOutputStream.cpp
 * Author: Mehmood Alam Khan
 * email: malagori@kth.se
 */

#include "BinaryDmOutputStream.hpp"
#include <cstdio>
#include <math.h>
#include <string>

using namespace std;

BinaryDmOutputStream::BinaryDmOutputStream(char *filename):DataOutputStream(filename) {
	if(filename != NULL) {
		writeToCout = false;
		delete fp;
		fp = NULL;
		file_was_opened = false;
		ofs = open_write_binary(filename);
	} else {
		ofs = &cout;
		writeToCout = true;
	}
}

BinaryDmOutputStream::~BinaryDmOutputStream() {
	if(ofs != NULL && !writeToCout) {
		delete ofs;
		ofs = NULL;
	}
}

void BinaryDmOutputStream::print(StrDblMatrix &dm) {
	float f;

	const size_t numNodes = dm.getSize();
	for ( size_t i=0; i<numNodes ; i++)
		cerr<<"row";
		for ( size_t j=i; j<numNodes; j++) {
			f=dm.getDistance(i,j);
			if (!isfinite(f))
				f=-1.0;
			cerr<<"\t"<<f<<endl;
			ofs->write(reinterpret_cast<char*>(&f),sizeof(f));
			}
	}

void BinaryDmOutputStream::printHeader( size_t numNodes ) {
	//converter variable is needed for running the binary output/input
	//also on 64-bit systems

	string tag = "FASTPHYLO 1";
	ofs->write(tag.c_str(), tag.length());
	long converter = numNodes;
	ofs->write(reinterpret_cast<char*>(&converter),sizeof(converter));
	//int nameSize=10;
	for(int i = 0; i < numNodes; ++i) {
		string name = m_names[i];
		cerr<<"name:\t"<<name<<endl;
		ofs->write(name.c_str(), name.length());
		ofs->write(":",1);
	}
}

void BinaryDmOutputStream::printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos) {
	m_names = names;
}
