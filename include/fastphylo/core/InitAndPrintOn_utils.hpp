//--------------------------------------------------
//
// File: InitAndPrintOn_utils.hpp
//
// Author: Isaac Elias
// e-mail: isaac@nada.kth.se
//
// cvs: $Id: InitAndPrintOn_utils.hpp,v 1.10 2006/12/19 09:24:36 isaac Exp $
//
//--------------------------------------------------
#pragma once

#include <string>
#include <iostream>
#include "fastphylo/core/file_utils.hpp"
#include "fastphylo/core/stl_utils.hpp"
#include "fastphylo/core/Sequence.hpp"

//
// This file contains function objects for printing and reading
// some common data types from streams.
//

//---------------------------------
// PRIMITIVE DATA TYPES
//
// Works for all data types for which the "<< i" and "i <<" operators
// are defined.
//
// Data_init<int>
// Data_init<float>
// Data_init<double>
// Data_init<std::string>
// Data_init<Sequence_float>

// INIT
template <class DataType> struct Data_init
{
    void operator()(std::istream &in, DataType &d) const
    {
        skipWhiteSpace(in);
        in >> d;
    }
};
template <class DataType> struct empty_Data_init
{
    void operator()(std::istream &in, const DataType &d) const
    {
        PROG_ERROR("Empty init function called");
    }
};

// PRINT ON
template <class DataType> struct Data_printOn
{

    std::ostream &operator()(std::ostream &os, const DataType &i) const
    {
        os << i;
        return os;
    }
};
template <class DataType> struct empty_Data_printOn
{

    std::ostream &operator()(std::ostream &os, const DataType &i) const
    {
        return os;
    }
};

//--------------------------------------
// SEQUENCE DOUBLE PAIR
//
struct Sequence_double
{
    Sequence s;
    double dbl;
};

std::istream &operator>>(std::istream &in, Sequence_double &strflt);

std::ostream &operator<<(std::ostream &os, const Sequence_double &strflt);
