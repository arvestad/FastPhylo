//--------------------------------------------------
//
// File: InitAndPrintOn_utils.cpp
//
// Author: Isaac Elias
// e-mail: isaac@nada.kth.se
//
// cvs: $Id: InitAndPrintOn_utils.cpp,v 1.6 2006/12/08 11:09:12 isaac Exp $
//
//--------------------------------------------------

#include "fastphylo/core/InitAndPrintOn_utils.hpp"
#include "fastphylo/core/Exception.hpp"
#include <cstdlib>
#include <string>
#include <cstring>
#include "fastphylo/core/file_utils.hpp"
#include "fastphylo/core/xml_output_global.hpp"

namespace
{

// atof()/atoi() replacements that, unlike their C-library namesakes,
// report a genuinely non-numeric token instead of silently treating it
// as 0. An empty token (a legitimate ":," / ":;" with nothing between,
// seen throughout this file's "str:flt"/"str:int" parsing) still parses
// as 0, matching the pre-existing atof()/atoi() behavior on "".
double parse_double_or_throw(const std::string &s)
{
    if (s.empty())
    {
        return 0;
    }
    char *endptr = nullptr;
    double val = std::strtod(s.c_str(), &endptr);
    if (endptr == s.c_str())
        THROW_EXCEPTION("expected a floating-point value, got \"" << s << "\"");
    return val;
}

int parse_int_or_throw(const std::string &s)
{
    if (s.empty())
    {
        return 0;
    }
    char *endptr = nullptr;
    long val = std::strtol(s.c_str(), &endptr, 10);
    if (endptr == s.c_str())
        THROW_EXCEPTION("expected an integer value, got \"" << s << "\"");
    return static_cast<int>(val);
}

} // namespace

std::istream &operator>>(std::istream &in, Sequence_double &strflt)
{
    // There are three cases:  1. "str:flt" 2. "str" 3. ":flt"
    skipWhiteSpace(in);

    // get the first part
    strflt.s = Sequence();
    while (strchr(":),;", in.peek()) == nullptr && in.peek() != EOF)
    {
        strflt.s.name += static_cast<char>(in.get());
    }

    // get the second part
    if (in.peek() == ':')
    {
        std::string fltstr;
        in.get(); // skip :
        while (strchr("), ;", in.peek()) == nullptr && in.peek() != EOF)
        {
            fltstr += static_cast<char>(in.get());
        }
        strflt.dbl = parse_double_or_throw(fltstr);
    }
    else
    {
        strflt.dbl = -1;
    }

    // if the last semicol of the tree
    if (in.peek() == ';')
    {
        in.get();
    }
    return in;
}

std::ostream &operator<<(std::ostream &os, const Sequence_double &strflt)
{
    if (xmlPrint)
    {
        if (strflt.dbl != -1)
        {
            os << " length=\"" << strflt.dbl << "\"";
        };
        os << ">";
        strflt.s.printShort(os);
    }
    else
    {
        strflt.s.printShort(os);
        if (strflt.dbl != -1)
        {
            os << ":" << strflt.dbl;
        }
    }
    return os;
}

std::istream &operator>>(std::istream &in, string_double &strflt)
{
    // There are three cases:  1. "str:flt" 2. "str" 3. ":flt"
    skipWhiteSpace(in);

    // get the first part
    strflt.s = std::string();
    while (strchr(":),;", in.peek()) == nullptr && in.peek() != EOF)
    {
        strflt.s += static_cast<char>(in.get());
    }

    // get the second part
    if (in.peek() == ':')
    {
        std::string fltstr;
        in.get(); // skip :
        while (strchr("), ;", in.peek()) == nullptr && in.peek() != EOF)
        {
            fltstr += static_cast<char>(in.get());
        }
        strflt.dbl = parse_double_or_throw(fltstr);
    }
    else
    {
        strflt.dbl = -1;
    }

    // if the last semicol of the tree
    if (in.peek() == ';')
    {
        in.get();
    }
    return in;
}

std::ostream &operator<<(std::ostream &os, const string_double &strflt)
{
    os << strflt.s;
    if (strflt.dbl != -1)
    {
        os << ":" << strflt.dbl;
    }
    return os;
}

std::istream &operator>>(std::istream &in, int_double &intdbl)
{
    // There are three cases:  1. "str:flt" 2. "str" 3. ":flt"
    skipWhiteSpace(in);

    // get the first part
    std::string s;
    while (strchr(":),;", in.peek()) == nullptr && in.peek() != EOF)
    {
        s += static_cast<char>(in.get());
    }
    intdbl.i = parse_int_or_throw(s);

    // get the second part
    if (in.peek() == ':')
    {
        std::string fltstr;
        in.get(); // skip :
        while (strchr("), ;", in.peek()) == nullptr && in.peek() != EOF)
        {
            fltstr += static_cast<char>(in.get());
        }
        intdbl.dbl = parse_double_or_throw(fltstr);
    }
    else
    {
        intdbl.dbl = -1;
    }

    // if the last semicol of the tree
    if (in.peek() == ';')
    {
        in.get();
    }
    return in;
}

std::ostream &operator<<(std::ostream &os, const int_double &intdbl)
{
    os << intdbl.i;
    if (intdbl.dbl != -1)
    {
        os << ":" << intdbl.dbl;
    }
    return os;
}

std::istream &operator>>(std::istream &in, string_int &strint)
{
    // There are three cases:  1. "str:flt" 2. "str" 3. ":flt"
    skipWhiteSpace(in);

    // get the first part
    strint.s = "";
    while (strchr(":),;", in.peek()) == nullptr && in.peek() != EOF)
    {
        strint.s += static_cast<char>(in.get());
    }

    // get the second part
    if (in.peek() == ':')
    {
        std::string fltstr;
        in.get(); // skip :
        while (strchr("), ;", in.peek()) == nullptr && in.peek() != EOF)
        {
            fltstr += static_cast<char>(in.get());
        }
        strint.i = parse_int_or_throw(fltstr);
    }
    else
    {
        strint.i = -1;
    }

    // if the last semicol of the tree
    if (in.peek() == ';')
    {
        in.get();
    }
    return in;
}

std::ostream &operator<<(std::ostream &os, const string_int &strint)
{
    os << strint.s;
    if (strint.i != -1)
    {
        os << ":" << strint.i;
    }
    return os;
}
