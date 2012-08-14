/*
 * BinaryDmOutputStream.cpp
 *
 *  Created on: Dec 1, 2011
 *      Author: Mehmood Alam Khan
 */

#include "BinaryDmOutputStream.hpp"

#include <cstdio>
#include <math.h>
#include <string>


using namespace std;

void
BinaryDmOutputStream::printRow( StrFloRow & dm, string name, int row) {
	int entriesPerRow = dm.getColumns();

	for( size_t j = row ; j < entriesPerRow ; j++ ) {
		float f = dm.getDistance(j);

		if ( ! isfinite(f) ){
			USER_WARNING("warning float not finite (use fix factor) " << f );
			f = -1.0;
		}

		ofs->write(reinterpret_cast<char*>( &f), sizeof f);
	}
}

void
BinaryDmOutputStream::printHeader( size_t numNodes ) {

	//dala mung da num tapa pa  lagawoo

	//converter variable is needed for running the binary output/input
	//also on 64-bit systems
	std::string tag = "FASTPHYLO 1";
	ofs->write(tag.c_str(), tag.length());
	long converter = numNodes;
	ofs->write(reinterpret_cast<char*>( &converter ), sizeof converter);
	//int nameSize=10;
	for(int i = 0; i < numNodes; ++i) {
		std::string name = m_names[i];
		ofs->write(name.c_str(), name.length());
		ofs->write("0",1);
	}

}

void
BinaryDmOutputStream::printStartRun(std::vector<std::string> & names, std::string & runId, Extrainfos &extrainfos) {
	m_names = names;
}
