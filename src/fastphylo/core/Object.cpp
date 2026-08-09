///////////////////////////////////////////////
//
// File: Object.cpp
//
// Author: Isaac Elias
//
///////////////////////////////////////////////

#include "fastphylo/core/Object.hpp"
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

/*===============================================
 * Object methods
 */

std::string Object::toString() const
{
    std::ostringstream out;
    printOn(out);
    return out.str();
}

ostream &Object::printOn(ostream &os) const
{
    // return os << "Object " << ((ptrdiff_t) this)<<" "; // if problem, then uncoment this line and comment below one
    return os << "Object " << (this) << " ";
}

/*=================================================
 * Out/In stream operators
 */
ostream &operator<<(ostream &os, const Object &obj)
{
    return obj.printOn(os);
}

ostream &operator<<(ostream &os, const Object *obj)
{
    return obj->printOn(os);
}

std::istream &operator>>(std::istream &in, Object *obj)
{
    return obj->objInitFromStream(in);
}
