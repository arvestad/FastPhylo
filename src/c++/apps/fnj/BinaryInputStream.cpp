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
// Templated over the matrix type: StrFloMatrix and StrDblMatrix expose
// the same setIdentifier()/getIdentifier()/getSize() interface, and the
// on-disk identifiers table doesn't depend on the distance value's type.
template <class Matrix>
void readIdentifiers(istream &fp, int newSize, Matrix &dm, std::vector<std::string> &names) {
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
// BinaryDmOutputStream::print()/printRow() always write raw `float`s on
// disk regardless of the matrix's own DistanceType (double for
// StrDblMatrix, float for StrFloMatrix) - setDistance() converts.
template <class Matrix>
bool readDistanceValues(istream &fp, int newSize, Matrix &dm) {
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

void BinaryInputStream::readHeaderIfNeeded() {
	if (input_was_read) {
		return;
	}
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

readstatus BinaryInputStream::readDM(StrFloMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos) {
	readHeaderIfNeeded();
	dm.resize(newSize);
	if (!input_was_read) {
		readIdentifiers(*fp, newSize, dm, names);
	}
	else {
		for (size_t i=0; i<names.size(); i++) {
			dm.setIdentifier(i,names.at(i));
}
		}
	if (!readDistanceValues(*fp, newSize, dm)) {
		return END_OF_RUN;
	}
	input_was_read=true;
	return DM_READ;
}

// This was an unimplemented stub ("Not implemented!", exit(-1)) - see
// fnj_binary_input_gap. This is the overload fnj's main.cpp actually
// calls (DataInputStream's virtual interface only declares the
// StrDblMatrix form - see this class's own DataInputStream.hpp's
// 2016-06-14 comment); the StrFloMatrix overload above has no live
// caller anywhere in the codebase, but is kept in sync since it reads
// the identical on-disk format. Mirrors BinaryDmOutputStream::print()/
// printHeader()'s write order exactly: one "FASTPHYLO 1"+size+names
// header per run, then one or more matrix bodies (original data, then
// any bootstrap replicates) until the stream runs out.
readstatus BinaryInputStream::readDM(StrDblMatrix &dm, vector<string> & names, string & runId, Extrainfos &extrainfos) {
	readHeaderIfNeeded();
	dm.resize(newSize);
	if (!input_was_read) {
		readIdentifiers(*fp, newSize, dm, names);
	}
	else {
		for (size_t i=0; i<names.size(); i++) {
			dm.setIdentifier(i,names.at(i));
		}
	}
	if (!readDistanceValues(*fp, newSize, dm)) {
		return END_OF_RUN;
	}
	input_was_read=true;
	return DM_READ;
}
