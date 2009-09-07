///////////////////////////////////////////////
//
// File: Object.cpp
//
// Author: Isaac Elias
//
// cvs: $Id: Object.cpp,v 1.6 2006/12/08 11:09:12 isaac Exp $
///////////////////////////////////////////////

#include "Object.hpp"
#include <string>
#include <stdio.h>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;


/*===============================================
 * Object methods
 */

std::string Object::toString() const{
  std::ostringstream out;
  printOn(out);
  return out.str();
}

ostream& Object::printOn(ostream& os) const{
  return os << "Object " << ((long) this)<<" ";
}


    
/*=================================================
 * Out/In stream operators
 */
ostream& operator<<(ostream & os, const Object *obj) {
    return obj->printOn(os);
}

ostream& operator<<(ostream & os, const Object &obj){
    return obj.printOn(os);
}

std::istream& operator>>(std::istream & in, Object &obj){
  return obj.objInitFromStream(in);
}

std::istream& operator>>(std::istream & in, Object *obj){
  return obj->objInitFromStream(in);
}


/*===================================================
 * String operators
 */
string operator + (const string & s, const Object & obj) {
    string t = s;
    t.append(obj.toString());
    return t;
}

string operator + (const string & s, const Object * obj) {
    string t = s;
    t.append(obj->toString());
    return t;
}

string operator + (const Object & obj, const string & s) {
    string t;
    t.append(obj.toString());
    t.append(s);
    return t;
}

string operator + (const Object * obj, const string & s) {
    string t;
    t.append(obj->toString());
    t.append(s);
    return t;
}




