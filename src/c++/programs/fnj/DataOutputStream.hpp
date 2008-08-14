#ifndef DATAOUTPUTSTREAM_HPP
#define DATAOUTPUTSTREAM_HPP


#include "Exception.hpp"
#include "SequenceTree.hpp"
#include "Extrainfos.hpp"

#include <iostream>
#include <vector>
#include <string>

class DataOutputStream
{
public:
  DataOutputStream( );
  DataOutputStream(char * filename );
  virtual ~DataOutputStream() {};
  // runId and extrainfos are just used by the XmlOutputStream
  virtual void print( tree2int_map & tree2count, bool noCounts, std::string & runId,  std::vector<std::string> & names, Extrainfos & extrainfos  ) = 0;
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
  virtual void print( tree2int_map & tree2count, bool noCounts, std::string & runId, std::vector<std::string> & names, Extrainfos & extrainfos );
};

#endif // DATAOUTPUTSTREAM_HPP
