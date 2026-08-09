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
#include "fastphylo/core/Object.hpp"
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

class Sequence : public Object
{
  public:
    std::string name;
    std::string seq;

    Sequence();
    Sequence(std::string n, std::string s);
    Sequence(const Sequence &s);

    Sequence &operator=(const Sequence &s);

    std::istream &objInitFromStream(std::istream &in) override;

    //------------------------------------------
    // PRINTING
    // Ex. "HUMAN    agct-agct"
    std::ostream &printOn(std::ostream &os) const override;
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

// The real formatting logic (Ex. "HUMAN    agct-agct"), found via ADL -
// object_modernization_plan.md Phase 1. printOn() above forwards here;
// this is what every ordinary `os << someSequence` call site actually
// binds to (an exact-match free function beats the inherited,
// conversion-requiring Object::operator<<).
std::ostream &operator<<(std::ostream &os, const Sequence &s);
