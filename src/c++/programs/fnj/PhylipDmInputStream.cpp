/*
 * PhylipDmInputStream.cpp
 *      Auther: Mehmood Alam Khan Email: malagori@kth.se
 */
#include "PhylipDmInputStream.hpp"
#include "DistanceMatrix.hpp"

#include <assert.h>

using namespace std;

PhylipDmInputStream::~PhylipDmInputStream() {
	if (file_was_opened)
		fin.close();
}

PhylipDmInputStream::PhylipDmInputStream(char *filename) {
	file_was_opened = false;
	if (filename ==NULL)
		fp = &cin;
	else {
		fin.open(filename, ifstream::in);
		if (fin.peek()==ifstream::traits_type::eof()) {
			fin.close();
			fin.clear();
			THROW_EXCEPTION("File is empty: \"" << filename << "\"");
		}
		if ( ! fin.good() ) {
			fin.close();
			fin.clear();
			THROW_EXCEPTION("File doesn't exist: \"" << filename << "\"");
		}
		file_was_opened = true;
		fp = & fin;
	}
}

//// check if file is empty or not
//bool isEmpty(std::ifstream& pFile)
//{
//	return pFile.peek() == std::ifstream::traits_type::eof();
//}

readstatus PhylipDmInputStream::readDM(StrDblMatrix &dm, vector<string> & names, string & runId, Extrainfos &extrainfos) {
	string line;
	int i1,i2,linePos,newSize;

	if (getline(*fp,line)==0)
		return END_OF_RUN;
	newSize=atoi(line.c_str());
	dm.resize(newSize);
	for (i1=0; i1<newSize; i1++) {
		getline(*fp,line);
		linePos=line.find_first_of(" \n\r\t");
		dm.setIdentifier(i1,line.substr(0,linePos));
		line=line.substr(line.find_first_not_of(" \n\r\t",linePos));
		for (i2=0; i2<newSize; i2++) {
			linePos=line.find_first_of(" \n\r\t");
			if (linePos>=0) {
				dm.setDistance(i1,i2,atof(line.substr(0,linePos).c_str()));
				line=line.substr(line.find_first_not_of(" \n\r\t",linePos));
				}
			else
				dm.setDistance(i1,i2,atof(line.c_str()));
			}
		}
	names.clear();
	for(size_t namei=0 ; namei<dm.getSize() ; namei++ )
		names.push_back(dm.getIdentifier(namei));
	return DM_READ;
}
