/*
 * BinaryDmOutputStream.cpp
 *
 *  Created on: March 18, 2013
 *      Author: Henric Zazzi
 */

#include "BinaryDmOutputStream.hpp"
#include <cstdio>
#include <math.h>
#include <string>
#include <vector>

using namespace std;

BinaryDmOutputStream::BinaryDmOutputStream(char *filename):DataOutputStream(filename) {
	if(filename != NULL) {
		writeToCout = false;
		// The base DataOutputStream(filename) constructor already opened fp
		// via fopen() (open_write_file()); we don't need it (we write
		// through ofs instead), so close it properly. This used to be
		// `delete fp` - deleting a FILE* that was never allocated with new -
		// which corrupted the heap and crashed on every invocation of
		// `-O binary -o <file>` (it only appeared to work when writing to
		// stdout, since that path never touches fp).
		fclose(fp);
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
	// One write() per row instead of one per float - same batching
	// rationale as PhylipDmOutputStream.cpp's printPHYLIPfastSD(), same
	// output speedup project (see phase0_audit.md). std::vector<float>
	// is contiguous/tightly-packed, so this writes the exact same bytes
	// in the exact same order as before, just in numNodes calls instead
	// of numNodes*(numNodes+1)/2.
	const size_t numNodes = dm.getSize();
	vector<float> row;
	row.reserve(numNodes);
	for ( size_t i=0; i<numNodes ; i++) {
		row.clear();
		for ( size_t j=i; j<numNodes; j++) {
			float f = dm.getDistance(i,j);
			if (!isfinite(f))
				f=-1.0;
			row.push_back(f);
		}
		ofs->write(reinterpret_cast<char*>(row.data()), row.size()*sizeof(float));
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
		ofs->write(name.c_str(), name.length());
		ofs->write(":",1);
	}
}

void BinaryDmOutputStream::printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos) {
	m_names = names;
}
