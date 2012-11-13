/*
 * PhylipDmInputStream.cpp
 *      Auther: Mehmood Alam Khan Email: malagori@kth.se
 */
#include "PhylipDmInputStream.hpp"
#include "DistanceMatrix.hpp"

using namespace std;

PhylipDmInputStream::~PhylipDmInputStream() {
	if ( file_was_opened )
		fin.close();
}

PhylipDmInputStream::PhylipDmInputStream(char * filename = 0 )  
{ 
	file_was_opened = false;
	if ( filename == 0 )
	{
		fp = & std::cin;    }
	else
	{
		fin.open(filename, ifstream::in );
		if (fin.peek() == std::ifstream::traits_type::eof())
		{
			fin.close();
			fin.clear();
			THROW_EXCEPTION("File is empty: \"" << filename << "\"");
		}

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
//// check if file is empty or not
//bool isEmpty(std::ifstream& pFile)
//{
//	return pFile.peek() == std::ifstream::traits_type::eof();
//}

readstatus
PhylipDmInputStream::readDM( StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos )
{
	dm.objInitFromStream(*fp);
	names.clear();

	for(size_t namei=0 ; namei<dm.getSize() ; namei++ ) {
		names.push_back(dm.getIdentifier(namei));
	}
	return DM_READ;
}

