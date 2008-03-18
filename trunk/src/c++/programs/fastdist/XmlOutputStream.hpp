#ifndef XMLOUTPUTSTREAM_HPP
#define XMLOUTPUTSTREAM_HPP

#include <cstdio>
#include "DataOutputStream.hpp"

class XmlOutputStream : public DataOutputStream
{
public:
  XmlOutputStream();
  XmlOutputStream(char * filename );
  virtual ~XmlOutputStream();
protected:
  virtual void print( StrDblMatrix & dm );
};

#endif // XMLOUTPUTSTREAM_HPP
