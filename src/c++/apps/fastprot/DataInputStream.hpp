#pragma once

#include <cstdio>
#include "fastphylo/core/Exception.hpp"
#include "fastphylo/core/Sequence.hpp"
#include "fastphylo/io/Extrainfos.hpp"

#include <iostream>
#include <fstream>

class DataInputStream
{
  public:
    DataInputStream()
    {
    }
    virtual ~DataInputStream()
    {
    }
    virtual bool read(std::vector<Sequence> &seqs, std::string &runId, std::vector<std::string> &names,
                      Extrainfos &extrainfos) = 0;
    virtual bool readSequences(std::vector<Sequence> &seqs, std::string &runId, Extrainfos &extrainfos) = 0;
};
