#include "FastaInputStream.hpp"
#include <cstdio>
#include <fstream>

using std::string;

FastaInputStream::~FastaInputStream() {
	if ( file_was_opened )
		fin.close();
}

FastaInputStream::FastaInputStream(char * filename = 0 )  
{ 
	file_was_opened = false;
	if ( filename == 0 )
	{
		fp = & std::cin;    }
	else
	{
		fin.open(filename, std::ifstream::in );
		if ( ! fin.good() )
		{
			fin.close();
			fin.clear();
			THROW_EXCEPTION("File doesn't exist: \"" << filename << "\"");
		}
		file_was_opened = true;
		fp = & fin;
	}
}



bool
FastaInputStream::read( std::vector<Sequence> &seqs, string & runId, std::vector<string> &names,  Extrainfos &extrainfos) 
{
	if ( ! readSequences(seqs,runId,extrainfos) ) return false;
	names.clear();
	names.reserve(seqs.size());
	for( size_t i=0;i<seqs.size();i++) {
		names.push_back(seqs[i].name);
	}
	return true;
}

bool
FastaInputStream::readSeq(std::vector<Sequence> &seqs, string &line, int linesRead) {

	linesRead++;
	seqs.resize(linesRead);
	string seqStr;
	Sequence &s = seqs[linesRead-1];
	string::size_type findPos = line.find_first_of("\t ",1);
	string::size_type nameEndPos;
	if ( findPos == string::npos ) {
		nameEndPos = line.size() -1;
	}
	else {
		nameEndPos = findPos-1;
	}

	s.name = line.substr(1, nameEndPos);

	bool readGreaterThan = false;
	while (  !readGreaterThan && getline ( *fp, line ) ) {
		//       c != 'n' && line[0] != '\0' && ! fp->eof() && fp->good()
		if ( line.size() > 0 ) {
			if ( line[0] == '>' ) {
				readGreaterThan = true;
				readSeq( seqs, line, linesRead );
			} else {
				seqStr+=line;
			}
		}
	}
// mehmood's changes here
	if (!fp->eof() && (seqStr.size() == 0 || seqStr.find_first_not_of("abcdefghiklmnopqrstuvwyzxABCDEFGHIKLMNOPQRSTUVWYZX -.?") != string::npos )) {
		THROW_EXCEPTION("Malformed Fasta format\n");
		exit(EXIT_FAILURE);
	}else {
		Sequence &s = seqs[linesRead-1];
		s.seq = seqStr;
	}
	return true;

}

bool
FastaInputStream::readSequences(std::vector<Sequence> &seqs,  std::string & runId, Extrainfos &extrainfos) {

	std::string line;
	while ( getline ( *fp, line ) ) {
		if ( line.size() > 0 && line[0] == '>' ) {
			readSeq( seqs , line ,0 );
			return true;
		}
	}
	THROW_EXCEPTION("Malformed Fasta format\n");
	exit(EXIT_FAILURE);
}
