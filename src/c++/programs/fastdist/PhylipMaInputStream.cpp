#include "PhylipMaInputStream.hpp"
#include <cstdio>

using namespace std;

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
PhylipMaInputStream::read( std::vector<std::string> &names, Extrainfos &extrainfos, std::vector<DNA_b128_String> &b128_strings) 
{
DNA_b128_StringsFromPHYLIP( *fp ,names,b128_strings);
 return true;
}

bool
PhylipMaInputStream::readSequences(std::vector<Sequence> &seqs, Extrainfos &extrainfos) {
  Sequence::readSequences(seqs ,*fp);
 return true;
}
