
#pragma once

#include <string>

#include <list>
// #include <ext/hash_map>
// #include <ext/hash_set>
#include <unordered_map>
#include <map>
#include <iostream>
#include <vector>
#include "fastphylo/core/file_utils.hpp"
#include "fastphylo/core/log_utils.hpp"
#include "fastphylo/core/Object.hpp"

//--------------------------
static inline void appendUntil(std::string &s, const char *chars, const int len, const char lastchar)
{
    int j = 0;
    for (; j < len; j++)
    {
        if (lastchar == chars[j])
            break;
    }
    s.append(chars, j);
}

static inline void appendAllNonChars(std::string &s, const char *chars, const int len, const char nonchar)
{
    int lastAppend = 0;
    for (int j = 0; j < len; j++)
    { // add all non space symbols (the last char is '\0')
        if (nonchar == chars[j])
        {
            s.append(chars + lastAppend, j - lastAppend);
            lastAppend = j + 1;
        }
    }
    s.append(chars + lastAppend, len - lastAppend);
}

//--------------------------
// vector printing

template <class T> std::istream &operator>>(std::istream &in, std::vector<T> &vec)
{

    vec.clear();

    skipWhiteSpace(in);

    if (in.peek() != '{')
    {
        USER_ERROR("Wrong vector format. e.g. {1,2,4}");
    }

    in.get(); // skip {

    while (in.peek() != '}')
    {
        T t;
        skipWhiteSpace(in);
        in >> t;
        vec.push_back(t);
        skipWhiteSpace(in);
        if (in.peek() != ',')
        {
            USER_ERROR("Wrong vector format. e.g. {1,2,4}");
        }
        in.get(); // skip ','
    }

    skipWhiteSpace(in);
    in.get(); // skip '}'
    return in;
}

template <class T> std::ostream &operator<<(std::ostream &os, std::vector<T> &vec)
{

    os << "{";
    for (size_t i = 0; i < vec.size() - 1; i++)
        os << vec[i] << " , ";

    os << vec[vec.size() - 1] << "}";
    return os;
}

//-----------------------------
// Hash table utils
struct eqstr
{
    bool operator()(const std::string &s1, const std::string &s2) const
    {
        return s1 == s2;
    }
};
struct hashstr
{
    std::hash<std::string> chhash;
    size_t operator()(const std::string &s) const
    {
        return chhash(s.c_str());
    }
};

// typedef __gnu_cxx::hash_map<const std::string, int, hashstr, eqstr> str2int_hashmap;
// typedef __gnu_cxx::hash_map<const std::string, std::string, hashstr, eqstr> str2str_hashmap;
using str2int_hashmap = std::unordered_map<const std::string, int, hashstr, eqstr>;
using str2str_hashmap = std::unordered_map<const std::string, std::string, hashstr, eqstr>;

//---------------------------------------------------------------
// Object 2 Object Map

using obj_ptr2obj_ptr_hashmap = std::unordered_map<Object *, Object *, objhash_ptr, objeq_ptr>;
using obj2obj_hashmap = std::unordered_map<Object, Object, objhash, objeq>;
