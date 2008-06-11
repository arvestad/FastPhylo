#include "DataOutputStream.hpp"
#include <cstdio>

using namespace std;

DataOutputStream::DataOutputStream(char * filename = 0)
{ 
  file_was_opened = false;
  if ( filename == 0 )
    {
   fp = & std::cout;    }
  else
    {
     fout.open(filename, ofstream::out );
     if ( ! fout.good() )
     {
       fout.close();
       fout.clear();
       THROW_EXCEPTION("File doesn't exist: \"" << filename << "\"");
     }
     file_was_opened = true;
     fp = & fout;
    }
}

TreeTextOutputStream::TreeTextOutputStream(char * filename = 0 ) : DataOutputStream( filename )
{ 
}


void
TreeTextOutputStream::print( tree2int_map & tree2count, bool noCounts ) 
{
  //OUTPUT THE TREES

  tree2int_map::iterator iter = tree2count.begin();
  for( ; iter!=tree2count.end() ; ++iter){
    if(noCounts)
      *fp << (*iter).first << endl;
    else
      *fp << (*iter).second << "  " << (*iter).first << endl; 
  }
}

