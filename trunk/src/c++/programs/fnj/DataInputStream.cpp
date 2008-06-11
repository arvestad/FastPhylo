#include "DataInputStream.hpp"
#include <cstdio>

using namespace std;

DataInputStream & chooseStream(char ** argv);

DataInputStream::DataInputStream()
{ 
}


PhylipMaInputStream::~PhylipMaInputStream() {
  if ( file_was_opened ) 
    fin.close();
}

PhylipMaInputStream::PhylipMaInputStream(char * filename = 0 )  
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

bool
PhylipMaInputStream::read( StrDblMatrix & dm, str2int_hashmap & name2id ) 
{
  dm.objInitFromStream(*fp);

  for(size_t namei=0 ; namei<dm.getSize() ; namei++ )
    {
      name2id[dm.getIdentifier(namei)] = namei;
    }


 return true;
}
