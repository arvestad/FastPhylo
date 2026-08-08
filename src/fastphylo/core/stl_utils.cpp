
#include "fastphylo/core/stl_utils.hpp"
#include <array>
#include <string>
#include <iostream>
#include "fastphylo/core/file_utils.hpp"

using namespace std;

string operator+(const string &s, const int i)
{
    // a 32bit int in decimal asciiz form is less than 12 characters
    std::array<char, 12> buf{};
#ifdef WIN32
    _snprintf(buf.data(), buf.size(), "%d", i);
#else
    snprintf(buf.data(), buf.size(), "%d", i);
#endif
    string t = s;
    t.append(buf.data());
    return t;
}

string operator+(const int i, const string &s)
{
    // a 32bit int in decimal asciiz form is less than 12 characters
    std::array<char, 12> buf{};
#ifdef WIN32
    _snprintf(buf.data(), buf.size(), "%d", i);
#else
    snprintf(buf.data(), buf.size(), "%d", i);
#endif
    string t = s;
    t.insert(0, buf.data());
    return t;
}

string operator+(const string &s, const float f)
{
    // a 32bit int in decimal asciiz form is less than 12 characters
    std::array<char, 20> buf{};
#ifdef WIN32
    _snprintf(buf.data(), buf.size(), "%f", f);
#else
    snprintf(buf.data(), buf.size(), "%f", f);
#endif
    string t = s;
    t.append(buf.data());
    return t;
}

string operator+(const float f, const string &s)
{
    // a 32bit int in decimal asciiz form is less than 12 characters
    std::array<char, 20> buf{};
#ifdef WIN32
    _snprintf(buf.data(), buf.size(), "%f", f);
#else
    snprintf(buf.data(), buf.size(), "%f", f);
#endif
    string t = s;
    t.insert(0, buf.data());
    return t;
}

std::istream &operator>>(std::istream &in, std::string *str)
{
    in >> *str;
    return in;
}

std::ostream &operator<<(std::ostream &os, std::string *str)
{
    os << *str;
    return os;
}

void print_map(str2int_map &m)
{

    cout << "map_size = " << m.size() << endl;
    for (const auto &kv : m)
    {
        cout << kv.first << " -> " << kv.second << endl;
    }
}
