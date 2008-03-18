#include "DataOutputStream.hpp"
#include <cstdio>

using namespace std;

DataOutputStream::DataOutputStream(char * filename = 0)
{ 
  fp = NULL;

  file_was_opened = false;
  if ( filename == 0 )
    {
      fp = stdout;
    }
  else
    {
      fp = open_write_file(filename);
    }
}

void
PhylipDmOutputStream::print( StrDblMatrix & dm ) 
{
  printPHYLIPfast(dm,fp,false);
}

