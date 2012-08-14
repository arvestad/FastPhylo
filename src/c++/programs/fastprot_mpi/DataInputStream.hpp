#ifndef DATAINPUTSTREAMM_HPP
#define DATAINPUTSTREAMM_HPP

#include <cstdio>
#include "../../Exception.hpp"
#include "../../Sequence.hpp"
#include "Extrainfos.hpp"



#include <iostream>
#include <fstream>



class DataInputStream
{
public:
  DataInputStream() {};
  virtual ~DataInputStream() {};
  virtual bool read( std::vector<Sequence> &seqs, std::string & runId, std::vector<std::string> &names, Extrainfos &extrainfos ) = 0;
  virtual bool readSequences( std::vector<Sequence> &seqs, std::string & runId, Extrainfos &extrainfos ) = 0;
};

#endif // DATAINPUTSTREAMM_HPP
