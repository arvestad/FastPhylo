#include "DataInputStream.hpp"
#include <cstdio>

using namespace std;

DataInputStream & chooseStream(char ** argv);

DataInputStream::DataInputStream()
{ 
}

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

readstatus
PhylipDmInputStream::readDM( StrDblMatrix & dm, std::vector<std::string> & names, Extrainfos & extrainfos ) 
{
 dm.objInitFromStream(*fp);
 names.clear();

 for(size_t namei=0 ; namei<dm.getSize() ; namei++ ) {
   names.push_back(dm.getIdentifier(namei));
 }
 return DM_READ;
}

