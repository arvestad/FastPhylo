/*
 * BinaryInputStream.hpp
 *
 * 	Created on: Dec 14, 2011
 *      Auther: Mehmood Alam Khan
 *       Email: malagori@kth.se
 */
#include "BinaryInputStream.hpp"
#include <cmath>
#include <cstring>

using namespace std;

BinaryInputStream::~BinaryInputStream() {
	if (file_was_opened) {
		fin.close();
}
}

BinaryInputStream::BinaryInputStream(char * filename)  {
	input_was_read=false;
	file_was_opened = false;
	if (filename==nullptr) { {
		fp = &cin;
	} } else {
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

namespace {
// Reads the ':'-terminated identifiers table once per binary file
// (input_was_read guards this in the caller, same as before this was
// extracted) and populates both dm's identifiers and names from it.
void readIdentifiers(istream &fp, int newSize, StrFloMatrix &dm, std::vector<std::string> &names) {
	char c;
	string identifier;
	for (int i=0; i< newSize; i++) {
		while(true){
			fp.read(&c, sizeof(c));
			if(c == ':') { break;
}
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
	for(size_t namei=0 ; namei<dm.getSize() ; namei++ ) {
		names.push_back(dm.getIdentifier(namei));
}
}

// Reads the newSize*(newSize+1)/2 upper-triangular distance values.
// Returns false on a short/failed read (caller returns END_OF_RUN).
bool readDistanceValues(istream &fp, int newSize, StrFloMatrix &dm) {
	for(int i = 0; i < newSize; ++i) {
		for(int j = i; j < newSize; ++j) {
			float f;
			if (!fp.read( reinterpret_cast<char*>( &f ), sizeof(f))) {
				return false;
}
			dm.setDistance(i, j, f);
		}
	}
	return true;
}
} // namespace

readstatus BinaryInputStream::readDM(StrFloMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos) {
	if (!input_was_read) {
		long converter;
		char t;
		int tagLength=11;
		string tag;
		for (int l=0; l< tagLength; ++l) {
			fp->read(&t, sizeof(t));
			tag+=t;
		}
		//converter variable is needed for running the binary output/input
		//also on 64-bit systems
		fp->read( reinterpret_cast<char*>( &converter ), sizeof(converter));
		newSize = static_cast<int>(converter);
	}
	dm.resize(newSize);
	if (!input_was_read) {
		readIdentifiers(*fp, newSize, dm, names);
	}
	else {
		for (int i=0; i<names.size(); i++) {
			dm.setIdentifier(i,names.at(i));
}
		}
	if (!readDistanceValues(*fp, newSize, dm)) {
		return END_OF_RUN;
	}
	input_was_read=true;
	return DM_READ;
}

readstatus BinaryInputStream::readDM(StrDblMatrix &dm, vector<string> & names, string & runId, Extrainfos &extrainfos) {
  std::cerr << "BinaryInputStream::readDM(StrDblMatrix, ...) -- Not implemented!" << endl;
  std::exit(-1);
}
