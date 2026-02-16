/*
 * BinaryInputStream.hpp
 *
 * 	Created on: Dec 14, 2011
 *      Auther: Mehmood Alam Khan
 *       Email: malagori@kth.se
 */
#include "BinaryInputStream.hpp"
#include <math.h>
#include <string.h>

using namespace std;

BinaryInputStream::~BinaryInputStream() {
	if (file_was_opened)
		fin.close();
}

BinaryInputStream::BinaryInputStream(char * filename)  {
	input_was_read=false;
	file_was_opened = false;
	if (filename==NULL)
		fp = &cin;
	else {
		fin.open(filename, ios::binary );
		if (!fin.good()) {
			fin.close();
			fin.clear();
			THROW_EXCEPTION("File doesn't exist: \"" << filename << "\"");
		}
		file_was_opened = true;
		fp = &fin;
	}
}

readstatus BinaryInputStream::readDM(StrFloMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos) {
	long converter;
	char t;
	int tagLength=11;
	string tag = "";

	if (!input_was_read) {
		for (int l=0; l< tagLength; ++l) {
			fp->read(&t, sizeof(t));
			tag+=t;
		}
		//converter variable is needed for running the binary output/input
		//also on 64-bit systems
		fp->read( reinterpret_cast<char*>( &converter ), sizeof(converter));
		newSize = (int)converter;
	}
	dm.resize(newSize);
	if (!input_was_read) {
		char c;
		string identifier = "";
		for (int i=0; i< newSize; i++) {
			while(true){
				fp->read(&c, sizeof(c));
				if(c == ':') break;
				identifier+=c;
			}
			// if identifier is empty don't add it to the sequence name vector
			if(identifier.empty()){
				continue;
			}
			dm.setIdentifier(i, identifier);
			identifier= "";
		}
		names.clear();
		for(size_t namei=0 ; namei<dm.getSize() ; namei++ )
			names.push_back(dm.getIdentifier(namei));
	}
	else {
		for (int i=0; i<names.size(); i++)
			dm.setIdentifier(i,names.at(i));
		}
	// read each line of the matrix and set the distances
	for(int i = 0; i < newSize; ++i) {
		for(int j = i; j < newSize; ++j) {
			float f;
			if (!fp->read( reinterpret_cast<char*>( &f ), sizeof(f)))
				return END_OF_RUN;
			dm.setDistance(i, j, f);
		}
	}
	input_was_read=true;
	return DM_READ;
}

readstatus BinaryInputStream::readDM(StrDblMatrix &dm, vector<string> & names, string & runId, Extrainfos &extrainfos) {
  std::cerr << "BinaryInputStream::readDM(StrDblMatrix, ...) -- Not implemented!" << endl;
  std::exit(-1);
}
