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
PhylipMaInputStream::read( std::vector<Sequence> &seqs, std::string & runId, std::vector<std::string> &names, Extrainfos &extrainfos )  
{
  //fungerar det här nedan verkligen?
  Sequence::readSequences(seqs, *fp);
  names.clear();names.reserve(seqs.size());
  for( size_t i=0;i<seqs.size();i++) {
    names.push_back(seqs[i].name);
  }
 return true;
}

bool
PhylipMaInputStream::readSequences(std::vector<Sequence> &seqs, std::string & runId, Extrainfos &extrainfos) {
  Sequence::readSequences(seqs ,*fp);
 return true;
}
