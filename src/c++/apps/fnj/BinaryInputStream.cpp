/*
 * BinaryInputStream.hpp
 *
 * 	Created on: Dec 14, 2011
 *      Auther: Mehmood Alam Khan
 *       Email: malagori@kth.se
 */
#include "BinaryInputStream.hpp"
#include "fastphylo/core/Exception.hpp"
#include <cmath>
#include <cstring>

using namespace std;

BinaryInputStream::~BinaryInputStream() {
	if (file_was_opened) {
		fin.close();
}
}

BinaryInputStream::BinaryInputStream(char * filename)  {
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
// Reads the ':'-terminated identifiers table (once per run - callers
// only invoke this right after a fresh header) and populates both dm's
// identifiers and names from it.
void readIdentifiers(istream &fp, int newSize, StrDblMatrix &dm, std::vector<std::string> &names) {
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
// disk regardless of the matrix's own DistanceType - setDistance()
// converts to double.
bool readDistanceValues(istream &fp, int newSize, StrDblMatrix &dm) {
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

bool BinaryInputStream::readNextHeader() {
	const int tagLength = 11;
	string tag;
	char t;
	for (int l = 0; l < tagLength; ++l) {
		if (!fp->read(&t, sizeof(t))) {
			return false; // true end of file - no more runs
		}
		tag += t;
	}
	// binary_dm_format_plan.md: the tag used to be read into a
	// variable and never actually checked - any garbage/foreign file
	// would silently misparse rather than fail clearly. "FASTPHYLO 2"
	// (bumped from "FASTPHYLO 1", same 11-byte length) is the version
	// that carries the matricesPerRun field read below.
	if (tag != "FASTPHYLO 2") {
		THROW_EXCEPTION("Bad binary distance-matrix tag: expected \"FASTPHYLO 2\", got \"" << tag << "\"");
	}
	//converter variable is needed for running the binary output/input
	//also on 64-bit systems
	long converter;
	fp->read( reinterpret_cast<char*>( &converter ), sizeof(converter));
	newSize = static_cast<int>(converter);
	long matricesPerRun;
	fp->read( reinterpret_cast<char*>( &matricesPerRun ), sizeof(matricesPerRun));
	bodiesRemainingInRun = static_cast<size_t>(matricesPerRun);
	return true;
}

// This was an unimplemented stub ("Not implemented!", exit(-1)) - see
// fnj_binary_input_gap. Mirrors BinaryDmOutputStream::printStartRun()/
// printHeader()/print()'s write order: one "FASTPHYLO 2"+size+
// matricesPerRun+names header per run, then exactly matricesPerRun
// matrix bodies (original data, then any bootstrap replicates) before
// the next run's header (or EOF) - binary_dm_format_plan.md's
// run-boundary fix. The sibling StrFloMatrix overload that used to
// exist here had zero callers anywhere in the codebase and was deleted
// rather than kept in sync with this state machine.
readstatus BinaryInputStream::readDM(StrDblMatrix &dm, vector<string> & names, string & runId, Extrainfos &extrainfos) {
	if (state == RunState::JustFinishedRun) {
		state = RunState::NeedHeader;
		return END_OF_RUN;
	}
	if (state == RunState::NeedHeader) {
		if (!readNextHeader()) {
			return END_OF_RUN; // true end of file
		}
	}
	dm.resize(newSize);
	if (state == RunState::NeedHeader) {
		readIdentifiers(*fp, newSize, dm, names);
		state = RunState::InRun;
	} else {
		for (size_t i=0; i<names.size(); i++) {
			dm.setIdentifier(i,names.at(i));
		}
	}
	if (!readDistanceValues(*fp, newSize, dm)) {
		// Truncated file: the header promised more bodies than actually
		// follow. Same "short read just means done" convention every
		// other reader here already uses, not a new failure mode.
		state = RunState::NeedHeader;
		return END_OF_RUN;
	}
	bodiesRemainingInRun--;
	if (bodiesRemainingInRun == 0) {
		state = RunState::JustFinishedRun;
	}
	return DM_READ;
}
