//--------------------------------------------------
//
// File: Sequence.hpp
//
// Author: Isaac Elias
// e-mail: isaac@nada.kth.se
//
//--------------------------------------------------
#pragma once

#include <string>
#include <vector>
#include "fastphylo/core/log_utils.hpp"
#include <fstream>
#include "fastphylo/core/file_utils.hpp"

//
// A simple class for handling named sequences.
//
// This class allows for reading and writing of PHYLIP sequence files
// and also bootstrapping.
//
//

class Sequence
{
  public:
    std::string name;
    std::string seq;

    Sequence();
    Sequence(std::string n, std::string s);
    Sequence(const Sequence &s);

    Sequence &operator=(const Sequence &s);

    //------------------------------------------
    // PRINTING
    // print only name
    virtual std::ostream &printShort(std::ostream &os) const;
    virtual void printWithoutGaps(std::ostream &os) const;

    //------------------------------------------
    // returns true if all characters in the sequences are in chars.
    bool onlyContains(std::string &chars);

    //----------------------------
    // READING PHYLIP SEQUENCEFILE
    // The vector is cleared of all sequences and new  sequences are added
    static void readSequences(std::vector<Sequence> &seqs, std::istream &in);
    static void printSequences(std::vector<Sequence> &seqs, std::ofstream &out);

    //-------------------------------------------
    // BOOTSTRAPPING
    // creates a bootstrapped data set from the sequences in seqs.
    static void bootstrapSequences(std::vector<Sequence> &seqs, std::vector<Sequence> &boot);
};

// Ex. "HUMAN    agct-agct". Found via ADL.
std::ostream &operator<<(std::ostream &os, const Sequence &s);
