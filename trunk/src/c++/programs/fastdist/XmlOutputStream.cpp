#include "XmlOutputStream.hpp"
#include <cstdio>
#include <libxml/xmlreader.h>

using namespace std;

XmlOutputStream::XmlOutputStream(char * filename = 0 ) : DataOutputStream(filename) 
{
  LIBXML_TEST_VERSION
  fprintf(fp,"<?xml version=\"1.0\"?>\n <root>\n  <runs>\n");
};

XmlOutputStream::~XmlOutputStream() 
{
  fprintf(fp,"  </runs>\n <root>\n");
};

void
XmlOutputStream::print( StrDblMatrix & dm ) 
{
  printPHYLIPfast(dm, fp , true );
}

