#ifndef DATAOUTPUTSTREAM_HPP
#define DATAOUTPUTSTREAM_HPP


#include "Exception.hpp"
#include "SequenceTree.hpp"

#include <iostream>

class DataOutputStream
{
public:
  DataOutputStream( );
  DataOutputStream(char * filename );
  virtual ~DataOutputStream() {};
  virtual void print( tree2int_map & tree2count, bool noCounts  ) = 0;
protected:
  std::ostream * fp;
  std::ofstream fout;
  bool file_was_opened;
};

class TreeTextOutputStream : public DataOutputStream
{
public:
  TreeTextOutputStream(char * filename );
  virtual ~TreeTextOutputStream() {};
  virtual void print( tree2int_map & tree2count, bool noCounts );
};

#endif // DATAOUTPUTSTREAM_HPP
