/*
 * BinaryInputStream.hpp
 *
 * 	Created on: Dec 14, 2011
 *      Auther: Mehmood Alam Khan
 *       Email: malagori@kth.se
 */
#include "BinaryInputStream.hpp"
#include <math.h>

using namespace std;

BinaryInputStream::~BinaryInputStream() {
	if ( file_was_opened )
		fin.close();
}

BinaryInputStream::BinaryInputStream(char * filename = 0 )  
{ 
	file_was_opened = false;
	//TODO: caused by changes in fastdist pipe will not work any more => TMP file?!?
	if ( filename == 0 )
	{
		fp = &std::cin;    }
	else
	{
		fin.open(filename, ios::binary );
		if ( ! fin.good() )
		{
			fin.close();
			fin.clear();
			THROW_EXCEPTION("File doesn't exist: \"" << filename << "\"");
		}
		file_was_opened = true;
		fp = &fin;
	}
}


readstatus
BinaryInputStream::readFloatDM( StrFloMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos )
{

	char t;
	int tagLength=11;
	std::string tag = "";
	for (int l=0; l< tagLength; ++l){
		fp->read(&t, sizeof(t));
		tag+=t;
	}
	//converter variable is needed for running the binary output/input
	//also on 64-bit systems
	long converter;
	fp->read( reinterpret_cast<char*>( &converter ), sizeof converter);

	int newSize = (int)converter;

	dm.resize(newSize);


	char c;
	std::string identifier = "";
	//printf("newSize:%d ", newSize);
	for (int i=0; i< newSize; i++){
		while(true){
			fp->read(&c, sizeof(c));
			if(c == '0') break;
			identifier+=c;
			//printf("%dname=%c\n",j,c);
		}
		dm.setIdentifier(i, identifier);
		identifier= "";
	}
	// read each line of the matrix and set the distances
	for(int i = 0; i < newSize; ++i) {
		for(int j = i; j < newSize; ++j) {
			float f;
			fp->read( reinterpret_cast<char*>( &f ), sizeof f);
			//printf("float=%f\n", f);
			dm.setDistance(i, j, f);
		}
	}

	names.clear();

	for(size_t namei=0 ; namei<dm.getSize() ; namei++ ) {
		names.push_back(dm.getIdentifier(namei));
	}

	return DM_READ;
}
